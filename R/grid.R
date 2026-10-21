# classydms: a common physical acquisition grid.
#
# Every sample of every cohort is resampled onto ONE grid whose rows are
# fixed RT seconds and whose columns are fixed CV volts, sized exactly to
# the model input (H x W). Pixel (i, j) then means the same
# (seconds, volts) in every image, no padding is needed, and positions
# can be carried in physical units end to end: peak_params and catalogs
# record the grid, and generators / heads map seconds / volts to pixels
# through it (R/placement.R).
#
# Why: trimming CV per cohort by signal occupancy gave each cohort a
# different column window (72-74 columns starting at different volts),
# so the same compound sat in different columns per cohort -- a
# cohort-specific (and so disease-correlated) cue. Some acquisitions also
# record 100 instead of 101 CV columns over the same voltage span, so CV
# must be resampled, not just cropped.

#' Define a common physical grid for GC-DMS images
#'
#' @param time Regular numeric vector of RT positions (seconds), one per
#'   image row.
#' @param cv Regular numeric vector of CV positions (volts), one per image
#'   column, increasing.
#' @param tol Relative tolerance for the regular-spacing check.
#' @return An object of class `dms_grid`: a list with `time`, `cv`, `H`,
#'   `W`, `dt`, `dcv`.
#' @seealso [make_dms_grid()], [resample_to_grid()]
#' @export
dms_grid <- function(time, cv, tol = 1e-6) {
  check_axis <- function(x, nm) {
    if (!is.numeric(x) || length(x) < 2L || anyNA(x))
      stop("dms_grid: `", nm, "` must be a numeric vector of length >= 2.")
    d <- diff(x)
    if (any(d <= 0) || diff(range(d)) > tol * abs(mean(d)))
      stop("dms_grid: `", nm, "` must be increasing and regularly spaced.")
    mean(d)
  }
  dt  <- check_axis(time, "time")
  dcv <- check_axis(cv, "cv")
  structure(list(time = as.numeric(time), cv = as.numeric(cv),
                 H = length(time), W = length(cv), dt = dt, dcv = dcv),
            class = "dms_grid")
}

#' Build a dms_grid from start, step and size
#'
#' @param H,W Grid rows / columns. Choose them as the model input size
#'   (H a multiple of 32, W a multiple of 4 for [dms_encoder()]); rows
#'   past a sample's last RT are filled with zeros by
#'   [resample_to_grid()].
#' @param t0,dt First RT (seconds) and RT step.
#' @param cv0,dcv First CV (volts) and CV step. Use [choose_cv_window()]
#'   to pick `cv0` from pooled signal occupancy.
#' @return A `dms_grid`.
#' @export
make_dms_grid <- function(H, W, t0 = 0, dt = 1, cv0, dcv = 0.55) {
  dms_grid(time = t0 + dt * (seq_len(H) - 1),
           cv   = cv0 + dcv * (seq_len(W) - 1))
}

#' @export
print.dms_grid <- function(x, ...) {
  cat(sprintf("dms_grid: %d x %d  RT %.2f..%.2f s (step %.3g)  CV %.3f..%.3f V (step %.3g)\n",
              x$H, x$W, x$time[1], x$time[x$H], x$dt,
              x$cv[1], x$cv[x$W], x$dcv))
  invisible(x)
}

# Linear interpolation matrix A (length(target) x length(src)) so that
# Z %*% t(A) interpolates each row of Z from `src` onto `target`.
# Targets outside [min(src), max(src)] get 0 (no signal), matching
# resample_to_time_grid()'s yleft = yright = 0. `method = "constant"`
# uses the nearest source column instead.
.axis_interp_matrix <- function(src, target, method = c("linear", "constant"),
                                tol = 1e-8) {
  method <- match.arg(method)
  if (is.unsorted(src)) stop("source axis must be increasing.")
  n_s <- length(src); n_t <- length(target)
  A <- matrix(0, n_t, n_s)
  span <- diff(range(src))
  inside <- target >= src[1] - tol * span & target <= src[n_s] + tol * span
  idx <- which(inside)
  if (!length(idx)) return(A)
  tt <- pmin(pmax(target[idx], src[1]), src[n_s])
  lo <- findInterval(tt, src, all.inside = TRUE)
  w  <- (tt - src[lo]) / (src[lo + 1L] - src[lo])
  if (method == "constant") {
    up <- w >= 0.5
    A[cbind(idx, ifelse(up, lo + 1L, lo))] <- 1
  } else {
    A[cbind(idx, lo)]      <- 1 - w
    A[cbind(idx, lo + 1L)] <- A[cbind(idx, lo + 1L)] + w
  }
  A
}

#' Resample one sample onto a CV grid
#'
#' Interpolates each RT row of `$Z` from the sample's own `$cv` axis onto
#' `target_cv`. Columns outside the sample's CV span are filled with 0.
#' Needed because acquisitions can differ in CV sampling (e.g. 100 vs
#' 101 columns over the same -40..15 V span).
#'
#' @param sample A sample list with `$cv` and `$Z` (`ncol(Z) == length(cv)`).
#' @param target_cv Increasing numeric vector of CV positions (volts).
#' @param method `"linear"` (default) or `"constant"` (nearest column).
#' @return The sample with `$Z` (`nrow(Z) x length(target_cv)`) and `$cv`
#'   set to `target_cv`.
#' @export
resample_to_cv_grid <- function(sample, target_cv,
                                method = c("linear", "constant")) {
  method <- match.arg(method)
  if (is.null(sample$cv) || is.null(sample$Z))
    stop("resample_to_cv_grid: sample must have $cv and $Z.")
  Z <- as.matrix(sample$Z); storage.mode(Z) <- "numeric"
  if (ncol(Z) != length(sample$cv))
    stop("resample_to_cv_grid: ncol(Z) (", ncol(Z),
         ") does not match length(cv) (", length(sample$cv), ").")
  src <- as.numeric(sample$cv)
  ord <- order(src)
  A <- .axis_interp_matrix(src[ord], target_cv, method)
  sample$Z  <- Z[, ord, drop = FALSE] %*% t(A)
  sample$cv <- as.numeric(target_cv)
  for (alt in c("Z_raw", "Z_smooth", "Z_pretrim")) sample[[alt]] <- NULL
  sample
}

#' Resample one sample onto a dms_grid (RT and CV)
#'
#' Applies [resample_to_time_grid()] with `grid$time`, then
#' [resample_to_cv_grid()] with `grid$cv`. The result is exactly
#' `grid$H x grid$W`, so every sample shares one physical pixel grid and
#' needs no padding.
#'
#' @param sample A sample list with `$time`, `$cv`, `$Z`.
#' @param grid A [dms_grid()].
#' @param method Interpolation method for both axes.
#' @return The resampled sample (`$time = grid$time`, `$cv = grid$cv`).
#' @export
resample_to_grid <- function(sample, grid, method = c("linear", "constant")) {
  stopifnot(inherits(grid, "dms_grid"))
  method <- match.arg(method)
  sample <- resample_to_time_grid(sample, grid$time, method = method)
  resample_to_cv_grid(sample, grid$cv, method = method)
}

#' Choose a common CV window from pooled signal occupancy
#'
#' Pools per-column signal occupancy (fraction of RT rows above each
#' sample's robust dust threshold, as in [trim_bounds_from_occupancy()])
#' across all supplied samples -- typically every cohort, excluding any
#' holdout -- on a regular CV axis, and returns the contiguous window of
#' `W` columns holding the most occupancy. The choice uses no labels.
#'
#' @param samples List of samples with `$cv` and `$Z` (any cohorts; CV
#'   sampling may differ, each is interpolated onto the common axis).
#' @param W Window width in columns (the model input width).
#' @param dcv CV step (volts) of the common axis.
#' @param cv_anchor A voltage the common axis passes through (keeps
#'   columns on the instrument's native positions). Default: the most
#'   common first CV value among the samples.
#' @return A list with `cv` (the `W` window voltages; pass `cv[1]` as
#'   `cv0` to [make_dms_grid()]), `occupancy` (pooled occupancy on the
#'   full common axis), `axis` (that axis), and `captured` (fraction of
#'   pooled occupancy inside the window).
#' @export
choose_cv_window <- function(samples, W, dcv = 0.55, cv_anchor = NULL) {
  stopifnot(length(samples) > 0L, W >= 1L)
  lo <- max(vapply(samples, function(s) min(s$cv), numeric(1)))
  hi <- min(vapply(samples, function(s) max(s$cv), numeric(1)))
  if (is.null(cv_anchor)) {
    firsts <- vapply(samples, function(s) min(s$cv), numeric(1))
    tab <- table(round(firsts, 6))
    cv_anchor <- as.numeric(names(tab)[which.max(tab)])
  }
  k_lo <- ceiling((lo - cv_anchor) / dcv - 1e-8)
  k_hi <- floor((hi - cv_anchor) / dcv + 1e-8)
  axis <- cv_anchor + dcv * (k_lo:k_hi)
  if (length(axis) < W)
    stop("choose_cv_window: the CV span shared by all samples has only ",
         length(axis), " columns at step ", dcv, "; W = ", W, " requested.")
  occ <- numeric(length(axis))
  for (s in samples) {
    rs <- resample_to_cv_grid(s, axis)
    eps <- robust_eps(rs$Z)
    occ <- occ + colMeans(rs$Z > eps)
  }
  occ <- occ / length(samples)
  sums <- stats::filter(occ, rep(1, W), sides = 1)[W:length(axis)]
  best <- which(sums == max(sums))
  # ties: prefer the most central window
  start <- best[which.min(abs(best - (length(axis) - W + 2) / 2))]
  win <- start:(start + W - 1L)
  list(cv = axis[win], occupancy = occ, axis = axis,
       captured = sum(occ[win]) / sum(occ))
}

# Pixel (1-based, fractional) of physical positions on a grid.
.grid_px <- function(grid, seconds, volts) {
  cbind((seconds - grid$time[1]) / grid$dt + 1,
        (volts   - grid$cv[1])   / grid$dcv + 1)
}
