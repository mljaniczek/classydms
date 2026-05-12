#' Read one raw GC-DMS data file
#'
#' Expects a tab-delimited file with this layout:
#'   - Line 1: literal "Vc"
#'   - Line 2: tab + compensation-voltage (CV) values, tab-separated
#'   - Line 3: header row (typically "Time Stamp" followed by per-CV labels)
#'   - Lines 4+: each row begins with a retention time (in seconds),
#'     followed by intensity values for each CV
#'
#' Implemented with base R only (no readr / stringr dependency) so it
#' works regardless of which packages the user has attached.
#'
#' @param path File path to one GC-DMS file.
#' @return List with `path`, `time` (numeric vector), `cv` (numeric vector),
#'   and `Z` (numeric matrix of intensity, rows = RT, cols = CV).
#' @export
read_dms_file <- function(path) {
  lines <- readLines(path, n = 3, warn = FALSE)
  if (length(lines) < 3) stop("File too short: ", path)

  # CV row is line 2; split on tab and drop empty tokens (leading tab)
  cv_tokens <- strsplit(lines[2], "\t", fixed = TRUE)[[1]]
  cv_tokens <- cv_tokens[nzchar(cv_tokens)]
  cv <- as.numeric(cv_tokens)

  # Data starts after 3 lines; read as TSV with no header
  dat <- utils::read.table(
    file = path,
    sep = "\t",
    skip = 3,
    header = FALSE,
    fill = TRUE,
    stringsAsFactors = FALSE,
    comment.char = ""
  )

  time <- as.numeric(dat[[1]])
  Z <- as.matrix(dat[, -1, drop = FALSE])
  storage.mode(Z) <- "numeric"

  # Strip any empty-string column/row names that read.table may have produced
  dimnames(Z) <- NULL

  if (anyNA(time) || anyNA(cv)) stop("Parsed NA in time/cv for: ", path)
  if (ncol(Z) != length(cv)) {
    stop("Intensity columns (", ncol(Z), ") != CV length (",
         length(cv), ") in: ", path)
  }

  list(path = path, time = time, cv = cv, Z = Z)
}

#' Read all GC-DMS files in a directory and assemble a samples list
#'
#' Convenience wrapper that lists files matching `pattern` in
#' `disease_dir` and `control_dir`, reads each with [read_dms_file()]
#' (skipping any that error), assigns each sample a name based on the
#' file basename so downstream results can be traced back to a specific
#' subject, and returns a parallel pair `samples` and `y`.
#'
#' Default `pattern` selects the positive ion channel only (POS files),
#' which is what the classification pipeline trains on. The negative
#' channel is not used.
#'
#' Sample names are set to the file basename with extension stripped
#' (e.g., `"241115-S1994-S1-DMS1_POS"`). This means you can index a
#' sample by name later --- e.g., `samples[["241115-S1994-S1-DMS1_POS"]]`
#' --- and so can `y` if you set `names(y) <- names(samples)`. The
#' function explicitly verifies that no name is the empty string,
#' which would otherwise propagate through `lapply()` and cause
#' "attempt to use zero-length variable name" errors downstream.
#'
#' @param disease_dir,control_dir Paths to directories.
#' @param pattern Regex matching the data files. Default selects
#'   positive-channel files; change to `"NEG.*\\.(xls|txt|tsv|csv)$"`
#'   for negative-channel-only.
#' @param verbose Whether to print progress messages.
#' @return List with `samples` (named list of sample lists) and `y`
#'   (named integer vector with 2 = disease, 1 = control).
#' @export
load_dms_directory <- function(disease_dir, control_dir,
                                pattern = "POS.*\\.(xls|txt|tsv|csv)$",
                                verbose = TRUE) {
  d_files <- list.files(disease_dir, pattern = pattern, full.names = TRUE)
  c_files <- list.files(control_dir, pattern = pattern, full.names = TRUE)
  if (verbose) {
    message("Found ", length(d_files), " disease file(s) and ",
            length(c_files), " control file(s)",
            " (pattern: ", pattern, ")")
  }
  if (length(d_files) == 0 || length(c_files) == 0)
    stop("No files found in one or both directories.")

  read_or_null <- function(p) {
    tryCatch(read_dms_file(p), error = function(e) {
      if (verbose) message("  Skipping ", basename(p), ": ",
                           conditionMessage(e))
      NULL
    })
  }
  samples_d <- lapply(d_files, read_or_null)
  samples_c <- lapply(c_files, read_or_null)

  keep_d <- !vapply(samples_d, is.null, logical(1))
  keep_c <- !vapply(samples_c, is.null, logical(1))
  samples_d <- samples_d[keep_d]
  samples_c <- samples_c[keep_c]
  d_files_kept <- d_files[keep_d]
  c_files_kept <- c_files[keep_c]

  samples <- c(samples_d, samples_c)

  # Assign meaningful names (file basename without extension) so each
  # sample is traceable back to a specific subject. Explicitly guard
  # against empty-string names, which can sneak in from c() / lapply on
  # filtered lists and would otherwise cause "attempt to use zero-length
  # variable name" errors downstream.
  sample_names <- tools::file_path_sans_ext(
    basename(c(d_files_kept, c_files_kept)))
  # De-duplicate any accidental name collisions
  if (anyDuplicated(sample_names) > 0) {
    sample_names <- make.unique(sample_names)
  }
  # Replace any empty / NA name with a positional fallback
  bad <- !nzchar(sample_names) | is.na(sample_names)
  if (any(bad)) sample_names[bad] <- paste0("sample_", which(bad))
  names(samples) <- sample_names

  y <- c(rep(2L, length(samples_d)), rep(1L, length(samples_c)))
  names(y) <- sample_names

  stopifnot(all(nzchar(names(samples))),
            length(samples) == length(y))

  if (verbose) {
    message("Loaded ", length(samples), " samples (",
            length(samples_d), " disease, ", length(samples_c), " control)")
    message("Sample names: ", paste(utils::head(sample_names, 3),
                                     collapse = ", "),
            if (length(sample_names) > 3)
              paste0(", ... (+", length(sample_names) - 3, " more)") else "")
  }
  list(samples = samples, y = y)
}


# classydms: preprocessing pipeline
# Smoothing, baseline correction, dust thresholding, occupancy-based
# trimming, normalization, and padding for GC-DMS intensity matrices.

#' Asymmetric Least Squares (AsLS) baseline correction
#'
#' Iteratively estimates a smooth baseline that lies under the data using


# classydms: preprocessing pipeline
# Smoothing, baseline correction, dust thresholding, occupancy-based
# trimming, normalization, and padding for GC-DMS intensity matrices.

#' Asymmetric Least Squares (AsLS) baseline correction
#'
#' Iteratively estimates a smooth baseline that lies under the data using
#' asymmetric weights so peaks are not absorbed.
#'
#' @param y Numeric vector (one RT trace).
#' @param lambda Smoothness penalty (default 1e6).
#' @param p Asymmetry parameter (default 0.01).
#' @param niter Number of iterations (default 10).
#' @return List with elements `z` (baseline) and `corrected` (`y - z`).
#' @keywords internal
asls_baseline <- function(y, lambda = 1e6, p = 0.01, niter = 10) {
  y <- as.numeric(y)
  m <- length(y)
  if (m < 3) return(list(z = y, corrected = rep(0, m)))
  D <- diff(diag(m), differences = 2)
  D <- Matrix::Matrix(D, sparse = TRUE)
  w <- rep(1, m)
  z <- y
  for (it in seq_len(niter)) {
    W <- Matrix::Diagonal(x = w)
    Zsol <- Matrix::solve(W + lambda * Matrix::t(D) %*% D, W %*% y)
    z <- as.numeric(Zsol)
    w <- ifelse(y > z, p, 1 - p)
  }
  list(z = z, corrected = y - z)
}

#' Savitzky-Golay smoothing + AsLS baseline correction (per CV column)
#'
#' Smooths each CV column along RT with a Savitzky-Golay filter, then
#' subtracts an AsLS-estimated baseline.
#'
#' @param Z RT x CV intensity matrix.
#' @param sg_p Polynomial order for Savitzky-Golay (default 3).
#' @param sg_n Window length (default 21; rounded up to odd, capped at nrow).
#' @param asls_lambda,asls_p,asls_niter AsLS parameters.
#' @param clip_floor Optional floor for corrected values (default 0).
#' @return List with `Z_smooth` (smoothed only) and `Z_corrected` (smoothed +
#'   baseline-corrected).
#' @export
preprocess_matrix_rt <- function(Z,
                                 sg_p = 3, sg_n = 21,
                                 asls_lambda = 1e6, asls_p = 0.01, asls_niter = 10,
                                 clip_floor = 0) {
  Z <- as.matrix(Z)
  storage.mode(Z) <- "numeric"
  nr <- nrow(Z); nc <- ncol(Z)

  if (sg_n %% 2 == 0) sg_n <- sg_n + 1
  sg_n <- min(sg_n, if (nr %% 2 == 1) nr else nr - 1)
  if (sg_n < 5) stop("sg_n too small after adjustment. Need >= 5.")

  Z_smooth <- Z; Z_corr <- Z
  for (j in seq_len(nc)) {
    y <- Z[, j]
    y_s <- signal::sgolayfilt(y, p = sg_p, n = sg_n)
    bl <- asls_baseline(y_s, lambda = asls_lambda, p = asls_p, niter = asls_niter)
    y_c <- bl$corrected
    if (!is.null(clip_floor)) y_c <- pmax(y_c, clip_floor)
    Z_smooth[, j] <- y_s
    Z_corr[, j] <- y_c
  }
  list(Z_smooth = Z_smooth, Z_corrected = Z_corr)
}

#' Run preprocessing on one sample (list with `$Z`, optional `$time`, `$cv`)
#'
#' Wraps [preprocess_matrix_rt()] and attaches the result back to the sample.
#'
#' @param sample A list with at least an element `Z` (RT x CV matrix).
#' @param sg_p,sg_n Savitzky-Golay parameters.
#' @param asls_lambda,asls_p,asls_niter AsLS parameters.
#' @param clip_floor Optional floor (default 0).
#' @return The sample with `Z_raw`, `Z_smooth`, and updated `Z`.
#' @export
process_one_sample <- function(sample,
                                sg_p = 3, sg_n = 21,
                                asls_lambda = 1e6, asls_p = 0.01, asls_niter = 10,
                                clip_floor = 0) {
  out <- preprocess_matrix_rt(Z = sample$Z,
                               sg_p = sg_p, sg_n = sg_n,
                               asls_lambda = asls_lambda, asls_p = asls_p,
                               asls_niter = asls_niter, clip_floor = clip_floor)
  sample$Z_raw <- sample$Z
  sample$Z_smooth <- out$Z_smooth
  sample$Z <- out$Z_corrected
  sample
}

#' Dust-threshold a sample below a fixed intensity
#'
#' Zeros out pixels of `sample$Z` whose intensity is below `basement_thr`.
#' Saves the pre-thresholded matrix into `sample$Z_pretrim` so peak-width
#' estimation can later use the unthresholded flanks.
#'
#' @param sample Sample list.
#' @param basement_thr Threshold below which pixels are zeroed.
#' @return Updated sample.
#' @export
baseline_basement <- function(sample, basement_thr) {
  Zb <- sample$Z
  Zb[Zb < basement_thr] <- 0
  sample$Z_pretrim <- sample$Z
  sample$Z <- Zb
  sample
}

#' Compute occupancy-based trim bounds for one sample
#'
#' For each RT row and CV column, computes the fraction of pixels above
#' the dust threshold (the "occupancy"), smooths these curves, and returns
#' the bounding box of rows/columns whose smoothed occupancy exceeds a
#' threshold.
#'
#' @param Z Intensity matrix (post-dust threshold).
#' @param eps Dust threshold (if NULL, computed per-sample as 1% of Q_0.95).
#' @param thr_rt,thr_cv Occupancy thresholds for retaining rows / columns.
#' @param smooth_k_rt,smooth_k_cv Rolling-mean window for occupancy curves.
#' @param min_keep_rt,min_keep_cv Minimum number of rows/columns to keep.
#' @return List with `rt_start`, `rt_end`, `cv_start`, `cv_end`, etc.
#' @export
trim_bounds_from_occupancy <- function(Z,
                                        eps = NULL,
                                        thr_rt = 0.005,
                                        thr_cv = 0.01,
                                        smooth_k_rt = 31,
                                        smooth_k_cv = 9,
                                        min_keep_rt = 100,
                                        min_keep_cv = 20) {
  Z <- as.matrix(Z); storage.mode(Z) <- "numeric"
  R <- nrow(Z); C <- ncol(Z)
  if (is.null(eps)) {
    pos <- Z[Z > 0]
    eps <- if (length(pos) < 50) 0 else as.numeric(stats::quantile(pos, 0.95)) * 0.01
  }
  M <- (Z > eps)
  occ_rt <- rowMeans(M); occ_cv <- colMeans(M)
  if (smooth_k_rt %% 2 == 0) smooth_k_rt <- smooth_k_rt + 1
  if (smooth_k_cv %% 2 == 0) smooth_k_cv <- smooth_k_cv + 1
  smooth_k_rt <- min(smooth_k_rt, if (R %% 2 == 1) R else R - 1)
  smooth_k_cv <- min(smooth_k_cv, if (C %% 2 == 1) C else C - 1)
  occ_rt_s <- zoo::rollmean(occ_rt, k = smooth_k_rt, fill = "extend", align = "center")
  occ_cv_s <- zoo::rollmean(occ_cv, k = smooth_k_cv, fill = "extend", align = "center")

  keep_rt <- which(occ_rt_s >= thr_rt)
  keep_cv <- which(occ_cv_s >= thr_cv)
  if (length(keep_rt) == 0 || length(keep_cv) == 0) {
    return(list(rt_start = 1, rt_end = R, cv_start = 1, cv_end = C,
                eps = eps, occ_rt = occ_rt_s, occ_cv = occ_cv_s))
  }
  rt_start <- min(keep_rt); rt_end <- max(keep_rt)
  cv_start <- min(keep_cv); cv_end <- max(keep_cv)
  if ((rt_end - rt_start + 1) < min_keep_rt) {
    mid <- floor((rt_start + rt_end) / 2)
    half <- floor(min_keep_rt / 2)
    rt_start <- max(1, mid - half)
    rt_end   <- min(R, rt_start + min_keep_rt - 1)
  }
  if ((cv_end - cv_start + 1) < min_keep_cv) {
    mid <- floor((cv_start + cv_end) / 2)
    half <- floor(min_keep_cv / 2)
    cv_start <- max(1, mid - half)
    cv_end   <- min(C, cv_start + min_keep_cv - 1)
  }
  list(rt_start = rt_start, rt_end = rt_end,
       cv_start = cv_start, cv_end = cv_end,
       eps = eps, occ_rt = occ_rt_s, occ_cv = occ_cv_s)
}

#' Apply trim bounds to a sample (clamps to actual size if smaller)
#'
#' @param sample A sample list (output of [read_dms_file()] or
#'   [process_one_sample()]), with elements `Z`, `time`, `cv`, and
#'   optionally `Z_pretrim`.
#' @param bounds A bounding box from [trim_bounds_from_occupancy()] or
#'   [pooled_trim_bounds()], with `rt_start`, `rt_end`, `cv_start`,
#'   `cv_end`.
#' @return The sample with its `Z` (and `Z_pretrim`, if present) trimmed
#'   to `bounds`, along with the corresponding subset of the `time` and
#'   `cv` axes. Indices are clamped to the actual matrix dimensions so a
#'   sample smaller than the trim bounds is not an error.
#' @export
trim_sample <- function(sample, bounds) {
  nr <- nrow(sample$Z); nc <- ncol(sample$Z)
  r1 <- max(1L, min(bounds$rt_start, nr))
  r2 <- max(r1,  min(bounds$rt_end,   nr))
  c1 <- max(1L, min(bounds$cv_start, nc))
  c2 <- max(c1,  min(bounds$cv_end,   nc))
  sample$time <- sample$time[r1:r2]
  sample$cv   <- sample$cv[c1:c2]
  sample$Z    <- sample$Z[r1:r2, c1:c2, drop = FALSE]
  if (!is.null(sample$Z_pretrim)) {
    nr_pre <- nrow(sample$Z_pretrim); nc_pre <- ncol(sample$Z_pretrim)
    r1_pre <- max(1L, min(r1, nr_pre)); r2_pre <- max(r1_pre, min(r2, nr_pre))
    c1_pre <- max(1L, min(c1, nc_pre)); c2_pre <- max(c1_pre, min(c2, nc_pre))
    sample$Z_pretrim <- sample$Z_pretrim[r1_pre:r2_pre, c1_pre:c2_pre, drop = FALSE]
  }
  sample
}

#' Pool per-sample trim bounds to a single cohort-wide bounding box
#'
#' Uses robust quantiles (5th of starts, 95th of ends by default) so a
#' single outlier sample with unusually narrow or wide content does not
#' force the cohort-wide bounds to lose data common to nearly all
#' samples.
#'
#' @param bounds_list List of bounds (one per sample) from
#'   [trim_bounds_from_occupancy()].
#' @param R,C Maximum row / column counts across samples; the pooled
#'   bounds are clamped to lie within `[1, R]` and `[1, C]`.
#' @param q_lo,q_hi Quantiles for pooling start / end positions
#'   (default 0.05 and 0.95).
#' @return List with `rt_start`, `rt_end`, `cv_start`, `cv_end` giving
#'   the cohort-wide bounding box.
#' @export
pooled_trim_bounds <- function(bounds_list, R, C, q_lo = 0.05, q_hi = 0.95) {
  rt_starts <- purrr::map_dbl(bounds_list, ~.x$rt_start)
  rt_ends   <- purrr::map_dbl(bounds_list, ~.x$rt_end)
  cv_starts <- purrr::map_dbl(bounds_list, ~.x$cv_start)
  cv_ends   <- purrr::map_dbl(bounds_list, ~.x$cv_end)
  list(
    rt_start = max(1, floor(stats::quantile(rt_starts, q_lo))),
    rt_end   = min(R, ceiling(stats::quantile(rt_ends,   q_hi))),
    cv_start = max(1, floor(stats::quantile(cv_starts, q_lo))),
    cv_end   = min(C, ceiling(stats::quantile(cv_ends,   q_hi)))
  )
}

#' Compute padding targets across a list of matrices
#'
#' Returns target H and W (rounded up to the nearest multiple of
#' `multiple`) so that all matrices in `Z_list` can be padded to a
#' common size. The `multiple = 32` default ensures clean compatibility
#' with the strided convolutions in the encoder (total stride along
#' RT is 32).
#'
#' @param Z_list List of 2-D numeric matrices.
#' @param multiple Integer multiple to round target dims up to
#'   (default 32).
#' @return List with `H` and `W` (the target dims), plus `raw_max_H`,
#'   `raw_max_W` (un-rounded maxes) and per-axis range vectors for
#'   diagnostics.
#' @export
compute_pad_targets <- function(Z_list, multiple = 32L) {
  dims <- do.call(rbind, lapply(Z_list, function(Z) {
    if (is.data.frame(Z)) Z <- as.matrix(Z); dim(Z)
  }))
  max_H <- max(dims[, 1]); max_W <- max(dims[, 2])
  target_H <- as.integer(ceiling(max_H / multiple) * multiple)
  target_W <- as.integer(ceiling(max_W / multiple) * multiple)
  list(H = target_H, W = target_W,
       raw_max_H = max_H, raw_max_W = max_W,
       range_H = range(dims[, 1]), range_W = range(dims[, 2]))
}

#' Center-pad a matrix to target dimensions
#'
#' Pads `Z` with `pad_value` so the result has dimensions
#' `target_H x target_W`. The original matrix is placed in the
#' center of the padded canvas. Errors if `Z` is already larger than
#' the target along either axis.
#'
#' @param Z A 2-D numeric matrix.
#' @param target_H,target_W Target row and column counts.
#' @param pad_value Fill value for the padded margins (default 0).
#' @return A `target_H x target_W` matrix.
#' @export
pad_to_target <- function(Z, target_H, target_W, pad_value = 0) {
  if (is.data.frame(Z)) Z <- as.matrix(Z)
  storage.mode(Z) <- "numeric"
  h <- nrow(Z); w <- ncol(Z)
  if (h > target_H || w > target_W) {
    stop("Matrix (", h, "x", w, ") exceeds target (", target_H, "x", target_W, ")")
  }
  if (h == target_H && w == target_W) return(Z)
  out <- matrix(pad_value, nrow = target_H, ncol = target_W)
  r_off <- floor((target_H - h) / 2)
  c_off <- floor((target_W - w) / 2)
  out[(r_off + 1):(r_off + h), (c_off + 1):(c_off + w)] <- Z
  out
}
