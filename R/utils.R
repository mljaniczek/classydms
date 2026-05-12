# classydms: utility functions
# Low-level helpers used across preprocessing, synthetic generation, and
# classification. Most are not exported; users typically interact with the
# higher-level wrappers in preprocess.R, pretrain.R, classify.R, etc.

#' Compute a per-sample relative "dust" threshold
#'
#' Returns 1% of the q-th quantile of positive values, used to identify
#' near-zero pixels that are likely numerical noise rather than real signal.
#' Returns 0 for samples with fewer than 50 positive values (treated as junk).
#'
#' @param Z A numeric matrix.
#' @param q Quantile to base the threshold on (default 0.95).
#' @param frac Fraction of the quantile to use (default 0.01).
#' @return A scalar threshold value.
#' @export
robust_eps <- function(Z, q = 0.95, frac = 0.01) {
  zpos <- Z[Z > 0]
  if (length(zpos) < 50) return(0)
  as.numeric(stats::quantile(zpos, q)) * frac
}

#' Measure peak width at a fraction of peak height (FWHM-style)
#'
#' Scans outward from a peak center along rows and columns until the
#' intensity drops below `frac_height` times the center value, returning
#' the width in pixels along each axis. Bounded by `max_radius_*` to prevent
#' walking into adjacent peaks.
#'
#' @param Z Intensity matrix.
#' @param r0,c0 Peak center indices (1-based).
#' @param frac_height Fraction of peak height defining the width boundary
#'   (default 0.5, i.e., FWHM).
#' @param eps Floor below which intensities are treated as zero.
#' @param max_radius_rt,max_radius_cv Maximum walk distance per axis.
#' @return Named numeric vector with `rt` and `cv` widths in pixels.
#' @keywords internal
peak_width_at <- function(Z, r0, c0, frac_height = 0.5, eps = 0,
                          max_radius_rt = Inf, max_radius_cv = Inf) {
  R <- nrow(Z); C <- ncol(Z)
  h0 <- Z[r0, c0]
  if (!is.finite(h0) || h0 <= max(eps, 0)) return(c(rt = NA_real_, cv = NA_real_))

  thr <- frac_height * h0

  r_lo <- r0
  while (r_lo > 1 && (r0 - r_lo + 1) < max_radius_rt &&
         Z[r_lo - 1, c0] >= thr) r_lo <- r_lo - 1
  r_hi <- r0
  while (r_hi < R && (r_hi - r0 + 1) < max_radius_rt &&
         Z[r_hi + 1, c0] >= thr) r_hi <- r_hi + 1
  w_rt <- r_hi - r_lo + 1

  c_lo <- c0
  while (c_lo > 1 && (c0 - c_lo + 1) < max_radius_cv &&
         Z[r0, c_lo - 1] >= thr) c_lo <- c_lo - 1
  c_hi <- c0
  while (c_hi < C && (c_hi - c0 + 1) < max_radius_cv &&
         Z[r0, c_hi + 1] >= thr) c_hi <- c_hi + 1
  w_cv <- c_hi - c_lo + 1

  c(rt = w_rt, cv = w_cv)
}

#' Greedily select top-K peak centers with spatial separation
#'
#' Returns up to `top_k` local-maximum pixels, ranked by intensity, subject
#' to a minimum separation constraint to prevent double-counting overlapping
#' peaks.
#'
#' @param Z Intensity matrix.
#' @param top_k Maximum number of peaks to return.
#' @param eps Threshold below which pixels are not considered peak candidates.
#'   If `NULL`, uses [robust_eps()].
#' @param min_sep_rt,min_sep_cv Minimum pixel separation between peaks.
#' @return Integer matrix with two columns (`r`, `c`) giving peak indices.
#' @keywords internal
pick_peak_centers <- function(Z, top_k = 30, eps = NULL,
                               min_sep_rt = 8, min_sep_cv = 2) {
  R <- nrow(Z); C <- ncol(Z)
  if (is.null(eps)) eps <- robust_eps(Z)

  idx <- which(Z > eps)
  if (length(idx) == 0) {
    return(matrix(integer(0), ncol = 2,
                  dimnames = list(NULL, c("r","c"))))
  }
  ord <- idx[order(Z[idx], decreasing = TRUE)]
  centers <- matrix(integer(0), ncol = 2,
                    dimnames = list(NULL, c("r","c")))

  for (k in seq_along(ord)) {
    if (nrow(centers) >= top_k) break
    rc <- arrayInd(ord[k], .dim = c(R, C))
    r0 <- rc[1]; c0 <- rc[2]
    if (nrow(centers) == 0) {
      centers <- rbind(centers, c(r0, c0))
      next
    }
    dr <- abs(centers[,1] - r0)
    dc <- abs(centers[,2] - c0)
    if (all(dr >= min_sep_rt | dc >= min_sep_cv)) {
      centers <- rbind(centers, c(r0, c0))
    }
  }
  centers
}

#' Log-quantile normalize a single sample
#'
#' Applies `log(1 + x) / Q_q(log(1 + x))`, where `Q_q` is the q-th quantile
#' of the log-transformed image. Falls back to dividing by 1 if the quantile
#' is below `eps` or non-finite. Per-sample normalization makes intensity
#' comparisons robust to instrument/subject-level absolute level differences.
#'
#' @param Z A non-negative numeric matrix.
#' @param q Quantile used as the normalization denominator (default 0.95).
#' @param eps Numerical floor on the denominator (default 1e-8).
#' @return A normalized matrix of the same shape as `Z`.
#' @export
normalize_sample <- function(Z, q = 0.95, eps = 1e-8) {
  Z <- log1p(Z)
  scale <- as.numeric(stats::quantile(Z, q))
  if (!is.finite(scale) || scale < eps) scale <- 1
  Z / scale
}

# ---- torch dataset helpers ----------------------------------------------

normalize_1based_indices <- function(idx, N, name = "indices") {
  idx <- as.integer(idx)
  if (length(idx) == 0) stop(name, " is empty.")
  if (anyNA(idx)) stop(name, " contains NA.")
  if (min(idx) == 0L && max(idx) <= (N - 1L)) idx <- idx + 1L
  if (min(idx) < 1L || max(idx) > N) {
    stop(name, " out of bounds. Range = [", min(idx), ", ", max(idx),
         "], expected within [1, ", N, "].")
  }
  idx
}

#' Build a torch dataset from a list of intensity matrices and labels
#'
#' Coerces matrices to a common (H, W) layout, stacks into a 4-D tensor of
#' shape `(N, 1, H, W)`, and returns a torch dataset whose items are
#' `list(x, y)` with one-channel images and integer labels.
#'
#' @param Z_list List of 2-D numeric matrices, all of the same dimensions
#'   (or transposes that will be auto-corrected).
#' @param y_vec Integer or factor labels.
#' @return An instantiated torch dataset.
#' @keywords internal
build_tensor_dataset <- function(Z_list, y_vec) {
  stopifnot(length(Z_list) == length(y_vec))
  if (length(Z_list) == 0) stop("Z_list is empty.")
  Z0 <- Z_list[[1]]
  if (is.data.frame(Z0)) Z0 <- as.matrix(Z0)
  if (!is.matrix(Z0)) Z0 <- as.matrix(Z0)
  storage.mode(Z0) <- "numeric"
  H <- nrow(Z0); W <- ncol(Z0)

  ensure_matrix_hw <- function(Z, H, W) {
    if (is.data.frame(Z)) Z <- as.matrix(Z)
    if (!is.matrix(Z)) Z <- as.matrix(Z)
    storage.mode(Z) <- "numeric"
    d <- dim(Z)
    if (all(d == c(H, W))) return(Z)
    if (all(d == c(W, H))) return(t(Z))
    stop("Matrix has dim ", paste(d, collapse = "x"),
         " but expected ", H, "x", W, " (or transposed).")
  }

  Z_list2 <- lapply(Z_list, ensure_matrix_hw, H = H, W = W)
  N <- length(Z_list2)
  X <- array(0, dim = c(N, 1, H, W))
  storage.mode(X) <- "double"
  for (i in seq_len(N)) X[i, 1, , ] <- Z_list2[[i]]

  X <- torch::torch_tensor(X, dtype = torch::torch_float())
  y <- torch::torch_tensor(as.integer(y_vec), dtype = torch::torch_long())

  ds_class <- torch::dataset(
    name = "dms_dataset",
    initialize = function() {},
    .getitem = function(i) list(x = X[i,..], y = y[i]),
    .length = function() X$size()[1]
  )
  ds_class()
}

#' View a subset of a torch dataset
#' @keywords internal
subset_dataset_view <- function(ds, indices) {
  N <- ds$.length()
  indices <- normalize_1based_indices(indices, N, name = "subset indices")
  ds_class <- torch::dataset(
    name = "subset_dataset_view",
    initialize = function() {},
    .getitem = function(i) ds$.getitem(indices[i]),
    .length = function() length(indices)
  )
  ds_class()
}

#' Stratified K-fold partition of binary labels
#' @keywords internal
make_stratified_folds <- function(y, k = 5, seed = NULL) {
  y <- as.integer(y)
  if (!is.null(seed)) set.seed(seed)
  cls <- sort(unique(y))
  if (length(cls) != 2) stop("Expected binary y with 2 unique values.")
  idx_a <- which(y == cls[1]); idx_b <- which(y == cls[2])
  fold_a <- sample(rep(1:k, length.out = length(idx_a)))
  fold_b <- sample(rep(1:k, length.out = length(idx_b)))
  folds <- vector("list", k)
  for (j in 1:k) {
    test_idx  <- c(idx_a[fold_a == j], idx_b[fold_b == j])
    train_idx <- setdiff(seq_along(y), test_idx)
    folds[[j]] <- list(train = train_idx, test = test_idx)
  }
  folds
}

#' AUC from binary labels and predicted probabilities
#' @keywords internal
auc_from_probs <- function(y_true, p_hat) {
  roc_obj <- pROC::roc(y_true, p_hat, quiet = TRUE, direction = "<")
  as.numeric(pROC::auc(roc_obj))
}

#' Per-batch on-the-fly augmentation (RT/CV shifts + intensity scaling + noise)
#' @keywords internal
augment_batch <- function(x,
                          rt_shift = 5,
                          cv_shift = 2,
                          intensity_scale = c(0.9, 1.1),
                          noise_sd = 0.0,
                          device = x$device) {
  sz <- x$size()
  B <- sz[[1]]; H <- sz[[3]]; W <- sz[[4]]

  if (!is.null(intensity_scale)) {
    a <- intensity_scale[1]; b <- intensity_scale[2]
    s <- torch::torch_rand(c(B, 1, 1, 1), device = device) * (b - a) + a
    x <- x * s
  }
  if (noise_sd > 0) {
    x <- x + torch::torch_randn_like(x) * noise_sd
    x <- x$clamp(min = 0)
  }
  if (rt_shift > 0 || cv_shift > 0) {
    dr <- if (rt_shift > 0) as.integer(sample(-rt_shift:rt_shift, B, replace = TRUE)) else rep(0L, B)
    dc <- if (cv_shift > 0) as.integer(sample(-cv_shift:cv_shift, B, replace = TRUE)) else rep(0L, B)
    x_aug <- torch::torch_zeros_like(x)
    for (i in 1:B) {
      r <- dr[i]; c <- dc[i]
      r_src1 <- max(1L, 1L - r); r_src2 <- min(H, H - r)
      r_dst1 <- max(1L, 1L + r); r_dst2 <- min(H, H + r)
      c_src1 <- max(1L, 1L - c); c_src2 <- min(W, W - c)
      c_dst1 <- max(1L, 1L + c); c_dst2 <- min(W, W + c)
      if (r_src1 <= r_src2 && c_src1 <= c_src2) {
        x_aug[i, 1, r_dst1:r_dst2, c_dst1:c_dst2] <-
          x[i, 1, r_src1:r_src2, c_src1:c_src2]
      }
    }
    x <- x_aug
  }
  x
}

#' Residual block used by encoder and tiny ResNet variants
#'
#' Two 3x3 convs with batch norm and a skip connection; optionally
#' downsamples via a 1x1 conv when `stride != 1` or channel counts differ.
#'
#' @keywords internal
#' @export
res_block <- torch::nn_module(
  initialize = function(in_ch, out_ch, stride = 1) {
    self$conv1 <- torch::nn_conv2d(in_ch, out_ch, kernel_size = 3, stride = stride,
                                    padding = 1, bias = FALSE)
    self$bn1   <- torch::nn_batch_norm2d(out_ch)
    self$conv2 <- torch::nn_conv2d(out_ch, out_ch, kernel_size = 3, stride = 1,
                                    padding = 1, bias = FALSE)
    self$bn2   <- torch::nn_batch_norm2d(out_ch)
    if (stride != 1 || in_ch != out_ch) {
      self$down <- torch::nn_sequential(
        torch::nn_conv2d(in_ch, out_ch, kernel_size = 1, stride = stride, bias = FALSE),
        torch::nn_batch_norm2d(out_ch)
      )
    } else {
      self$down <- NULL
    }
  },
  forward = function(x) {
    identity <- x
    out <- self$conv1(x) |> self$bn1() |> torch::nnf_relu()
    out <- self$conv2(out) |> self$bn2()
    if (!is.null(self$down)) identity <- self$down(identity)
    torch::nnf_relu(out + identity)
  }
)
