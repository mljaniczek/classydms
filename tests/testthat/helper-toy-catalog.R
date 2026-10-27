# Toy peak_catalog for generator tests: `n` compounds on a fraction-mode
# canvas, each with a handful of valid (sigma_rt, sigma_cv, intensity)
# observations and the given prevalence.
make_toy_catalog <- function(prevalence = c(0.95, 0.9, 0.6, 0.4, 0.25,
                                            0.1, 0.05, 0.8, 0.5, 0.3),
                             n_obs = 4L) {
  n <- length(prevalence)
  rep_obs <- function(v) lapply(seq_len(n), function(i) rep(v, n_obs))
  structure(list(
    compounds = tibble::tibble(
      compound_id   = seq_len(n),
      rt_loc        = seq(0.08, 0.92, length.out = n),
      cv_loc        = rep(c(0.3, 0.7), length.out = n),
      prevalence    = prevalence,
      intensity_obs = lapply(seq_len(n), function(i) seq(1, 2, length.out = n_obs)),
      sigma_rt_obs  = rep_obs(3),
      sigma_cv_obs  = rep_obs(1.5)),
    noise = list(mean = 0, sd = 0.01),
    parameters = list(coord_mode = "fraction",
                      data_dims = cbind(rows = 120, cols = 24))),
    class = "peak_catalog")
}
