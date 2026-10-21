# Raw-like samples on an instrument CV axis (-40..15 V, `n_cv` columns)
# and an irregular RT axis, with Gaussian peaks at physical positions.
raw_samples <- function(n, peaks, n_cv = 101L, seed = 1L) {
  set.seed(seed)
  cv <- seq(-40, 15, length.out = n_cv)
  lapply(seq_len(n), function(s) {
    time <- cumsum(c(0, stats::runif(179, 0.9, 1.3)))
    Z <- matrix(abs(stats::rnorm(length(time) * n_cv, 0, 1e-3)),
                length(time), n_cv)
    for (k in seq_len(nrow(peaks)))
      Z <- Z + peaks$amp[k] *
        outer(exp(-0.5 * ((time - peaks$sec[k]) / 2)^2),
              exp(-0.5 * ((cv - peaks$volt[k]) / 0.45)^2))
    list(time = time, cv = cv, Z = Z)
  })
}
grid_peaks <- data.frame(sec = c(30, 80, 150), volt = -40 + 0.55 * c(30, 55, 80),
                         amp = c(1, 2, 3))

test_that("dms_grid validates and maps physical units to pixels", {
  g <- make_dms_grid(H = 224, W = 32, t0 = 0, dt = 1, cv0 = -30, dcv = 0.55)
  expect_s3_class(g, "dms_grid")
  expect_equal(c(g$H, g$W), c(224L, 32L))
  expect_equal(.grid_px(g, 10, -30 + 0.55 * 4), cbind(11, 5))
  expect_error(dms_grid(c(0, 1, 3), 1:3), "regularly spaced")
})

test_that("CV resampling is exact on matching axes and aligns 100-col runs", {
  s <- raw_samples(1L, grid_peaks)[[1]]
  same <- resample_to_cv_grid(s, s$cv)
  expect_equal(same$Z, s$Z, tolerance = 1e-12)
  s100 <- raw_samples(1L, grid_peaks, n_cv = 100L)[[1]]
  target <- -40 + 0.55 * (0:100)
  r101 <- resample_to_cv_grid(s, target)
  r100 <- resample_to_cv_grid(s100, target)
  row <- which.min(abs(s$time - 150))
  expect_equal(which.max(r101$Z[row, ]), which.max(r100$Z[row, ]))
  # outside the source span -> 0
  expect_equal(resample_to_cv_grid(s, c(-50, -45))$Z[1, ], c(0, 0))
})

test_that("choose_cv_window keeps the occupied band", {
  smp <- c(raw_samples(3L, grid_peaks), raw_samples(3L, grid_peaks, 100L, 2L))
  win <- choose_cv_window(smp, W = 64L)
  expect_length(win$cv, 64L)
  expect_true(all(grid_peaks$volt >= min(win$cv) & grid_peaks$volt <= max(win$cv)))
  expect_equal(diff(win$cv)[1], 0.55)
  expect_error(choose_cv_window(smp, W = 200L), "only")
})

test_that("grid pipeline places peaks by seconds/volts across cohorts", {
  a_raw <- raw_samples(5L, grid_peaks, 101L, seed = 1L)
  b_raw <- raw_samples(5L, grid_peaks, 100L, seed = 2L)
  win <- choose_cv_window(c(a_raw, b_raw), W = 64L)
  g <- make_dms_grid(H = 224, W = 64, t0 = 0, dt = 1, cv0 = win$cv[1])
  to_pp <- function(raw, tag) {
    Z <- lapply(raw, function(s) resample_to_grid(s, g)$Z)
    names(Z) <- paste0(tag, seq_along(Z))
    estimate_peak_params(Z, grid = g, top_k = 3, verbose = FALSE)
  }
  a <- to_pp(a_raw, "a"); b <- to_pp(b_raw, "b")
  expect_identical(a$grid, g)
  # the same compound sits in the same column in both cohorts
  for (k in seq_len(nrow(grid_peaks))) {
    near <- function(pp) round(stats::median(pp$cv_loc_raw[
      abs(pp$rt_seconds_raw - grid_peaks$sec[k]) < 3] * g$W))
    expect_equal(near(a), near(b))
  }

  ct <- build_universal_catalog(list(a = a, b = b), eps_rt = 3, eps_cv = 1,
                                min_support_frac = 0.1,
                                coord_mode = "physical")
  expect_identical(ct$parameters$grid, g)
  cm <- ct$compounds[order(ct$compounds$rt_loc), ]
  expect_equal(nrow(cm), nrow(grid_peaks))
  truth_px <- .grid_px(g, grid_peaks$sec, grid_peaks$volt)

  pos <- catalog_pixel_positions(ct, g$H, g$W)
  pos <- pos[match(cm$compound_id, pos$compound_id), ]
  expect_true(all(abs(pos$rt_px - truth_px[, 1]) <= 1))
  expect_true(all(abs(pos$cv_px - truth_px[, 2]) <= 0.5))

  res <- generate_one_synthetic_from_catalog(
    ct, g$H, g$W, add_noise = FALSE, size_jitter = 0,
    location_jitter_rt = 0, location_jitter_cv = 0,
    prevalence_override = stats::setNames(rep(1, nrow(cm)),
                                          as.character(cm$compound_id)))
  pl <- res$placements[match(cm$compound_id, res$placements$compound_id), ]
  expect_true(all(abs(pl$mu_rt - truth_px[, 1]) <= 1))
  expect_true(all(abs(pl$mu_cv - truth_px[, 2]) <= 0.5))
  # each placement is exactly one observation's seconds/volts on the grid
  for (i in seq_len(nrow(cm))) {
    obs <- .grid_px(g, cm$rt_seconds_obs[[i]], cm$cv_volts_obs[[i]])
    expect_true(any(abs(obs[, 1] - pl$mu_rt[i]) < 1e-9 &
                    abs(obs[, 2] - pl$mu_cv[i]) < 1e-9))
  }

  geo <- compound_head_geometry(ct, g$H, g$W)
  expect_equal(geo$cv_frac[match(cm$compound_id, geo$compound_id)] * g$W + 0.5,
               pos$cv_px)
  # larger canvas: centred; smaller canvas: refused
  pos2 <- catalog_pixel_positions(ct, g$H + 32, g$W + 8)
  expect_equal(pos2$cv_px, catalog_pixel_positions(ct, g$H, g$W)$cv_px + 4)
  expect_error(catalog_pixel_positions(ct, g$H - 32, g$W), "smaller")
  # peak_params generator path uses the grid too
  px <- .place_pp_peaks(a, seq_along(a$rt_seconds_raw), g$H, g$W)
  expect_equal(px, .grid_px(g, a$rt_seconds_raw, a$cv_volts_raw))
})

test_that("grid consistency is enforced", {
  g <- make_dms_grid(H = 224, W = 64, cv0 = -30)
  raw <- raw_samples(2L, grid_peaks)
  Z <- lapply(raw, function(s) resample_to_grid(s, g)$Z)
  names(Z) <- c("x", "y")
  expect_error(estimate_peak_params(list(x = Z[[1]][1:200, ]), grid = g),
               "resample_to_grid")
  pp <- estimate_peak_params(Z, grid = g, top_k = 3, verbose = FALSE)
  bad <- pp; bad$data_dims[1, 2] <- 63L
  expect_error(validate_peak_params(bad), "grid")
  g2 <- make_dms_grid(H = 224, W = 64, cv0 = -29.45)
  pp2 <- pp; pp2$grid <- g2
  expect_message(m <- merge_peak_params(list(a = pp, b = pp2)), "same dms_grid")
  expect_null(m$grid)
  expect_identical(merge_peak_params(list(a = pp, b = pp))$grid, g)
  expect_identical(subset_peak_params(pp, 2L)$grid, g)
})
