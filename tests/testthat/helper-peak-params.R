# Toy peak_params generator for catalog-integrity tests.
#
# Builds a peak_params object shaped like estimate_peak_params() output
# with physical axes: a 1 s/row RT grid (time = 0, 1, ..., R - 1) and a
# 0.5 V/column CV grid. Each compound has a fixed true (row, col), its
# own intensity and sigma_rt (so per-observation fields can be traced
# back to the compound they came from), and a per-sample firing
# probability. `private` compounds fire in every listed sample and
# nowhere else, which lets tests check that excluded samples contribute
# nothing to a catalog.
make_toy_peak_params <- function(n_samples, compounds, R = 300L, C = 40L,
                                 cv_min = -20, private = NULL,
                                 seed = 1L) {
  set.seed(seed)
  time_axis <- seq(0, by = 1, length.out = R)
  cv_axis   <- seq(cv_min, by = 0.5, length.out = C)
  rows <- cols <- sig_rt <- sig_cv <- inten <- sidx <- integer(0)
  inten <- numeric(0); sig_rt <- numeric(0); sig_cv <- numeric(0)
  add_peak <- function(cmp, s) {
    r <- cmp$row + sample(-1:1, 1L)
    rows  <<- c(rows, r)
    cols  <<- c(cols, cmp$col)
    inten <<- c(inten, cmp$intensity * exp(stats::rnorm(1L, 0, 0.05)))
    sig_rt <<- c(sig_rt, cmp$sigma_rt)
    sig_cv <<- c(sig_cv, 1)
    sidx  <<- c(sidx, s)
  }
  for (s in seq_len(n_samples)) {
    for (k in seq_len(nrow(compounds))) {
      if (stats::runif(1L) < compounds$prevalence[k]) add_peak(compounds[k, ], s)
    }
    if (!is.null(private) && s %in% private$samples) {
      for (k in seq_len(nrow(private$compounds)))
        add_peak(private$compounds[k, ], s)
    }
  }
  n_per <- tabulate(sidx, nbins = n_samples)
  list(
    n_peaks   = list(values = n_per, mean = mean(n_per), sd = stats::sd(n_per)),
    rt_loc    = list(values = rows / R, mean = mean(rows / R), sd = stats::sd(rows / R)),
    cv_loc    = list(values = cols / C, mean = mean(cols / C), sd = stats::sd(cols / C)),
    noise     = list(mean = 0.01, sd = 0.005),
    n_samples = as.integer(n_samples),
    n_peaks_detected = length(rows),
    sigma_rt_raw   = sig_rt,
    sigma_cv_raw   = sig_cv,
    intensity_raw  = inten,
    rt_loc_raw     = rows / R,
    cv_loc_raw     = cols / C,
    sample_idx_raw = sidx,
    sample_names   = sprintf("s%03d", seq_len(n_samples)),
    rt_seconds_raw = time_axis[rows],
    cv_volts_raw   = cv_axis[cols],
    coord_mode     = "physical",
    time_ranges = do.call(rbind, rep(list(c(min = min(time_axis),
                                            max = max(time_axis), n = R)),
                                     n_samples)),
    cv_ranges   = do.call(rbind, rep(list(c(min = min(cv_axis),
                                            max = max(cv_axis), n = C)),
                                     n_samples)),
    data_dims   = cbind(rows = rep(R, n_samples), cols = rep(C, n_samples))
  )
}

toy_compounds <- function(n = 12L, R = 300L, C = 40L, seed = 2L) {
  set.seed(seed)
  data.frame(
    row       = round(seq(20, R - 20, length.out = n)),
    col       = sample(5:(C - 5), n, replace = TRUE),
    intensity = round(exp(stats::rnorm(n, 0, 1)), 3),
    sigma_rt  = seq_len(n) + 0.5,
    prevalence = stats::runif(n, 0.4, 0.95)
  )
}

# Real-looking samples for placement tests: Gaussian peaks rendered into
# R x C matrices on a 1 s/row RT axis (time = row - 1) and a 0.5 V/col
# CV axis starting at `cv_start`. `peaks` has true `row`, `col`.
make_toy_samples <- function(n, peaks, R = 200L, C = 36L, cv_start = -10,
                             seed = 1L) {
  set.seed(seed)
  rr <- seq_len(R); cc <- seq_len(C)
  Z <- lapply(seq_len(n), function(s) {
    m <- matrix(abs(stats::rnorm(R * C, 0, 1e-3)), R, C)
    for (k in seq_len(nrow(peaks)))
      m <- m + peaks$amp[k] * outer(exp(-0.5 * ((rr - peaks$row[k]) / 1.5)^2),
                                    exp(-0.5 * ((cc - peaks$col[k]) / 0.8)^2))
    m
  })
  names(Z) <- sprintf("s%02d_%d", seq_len(n), C)
  list(Z = Z,
       time = rep(list(seq(0, by = 1, length.out = R)), n),
       cv   = rep(list(seq(cv_start, by = 0.5, length.out = C)), n))
}
