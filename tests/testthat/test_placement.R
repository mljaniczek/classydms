toy_peaks <- data.frame(row = c(20, 70, 120, 180), col = c(3, 12, 25, 34),
                        amp = c(1, 2, 1.5, 3))

pp_from_samples <- function(smp) {
  estimate_peak_params(smp$Z, time_axes = smp$time, cv_axes = smp$cv,
                       top_k = 4, verbose = FALSE)
}

test_that(".pad_frac_to_px matches pad_to_target's centre placement", {
  Z <- matrix(0, 200, 36); Z[120, 25] <- 1
  P <- pad_to_target(Z, 224, 48)
  hit <- which(P == 1, arr.ind = TRUE)
  expect_equal(.pad_frac_to_px(120 / 200, 200, 224), unname(hit[1, "row"]))
  expect_equal(.pad_frac_to_px(25 / 36, 36, 48), unname(hit[1, "col"]))
  # no padding possible -> scale (resize semantics); unknown dims -> scale
  expect_equal(.pad_frac_to_px(0.5, 300, 224), 112)
  expect_equal(.pad_frac_to_px(0.5, NA, 48), 24)
})

test_that("estimate_peak_params records data dims and catalogs carry them", {
  pp <- pp_from_samples(make_toy_samples(6L, toy_peaks))
  expect_equal(unname(pp$data_dims[1, ]), c(200L, 36L))
  expect_equal(nrow(pp$data_dims), 6L)
  expect_silent(validate_peak_params(pp))
  s <- subset_peak_params(pp, c(2L, 5L))
  expect_equal(nrow(s$data_dims), 2L)
  ct <- build_peak_catalog(pp, eps_rt = 3, eps_cv = 1, coord_mode = "physical")
  expect_equal(ct$parameters$data_dims, pp$data_dims)
  expect_equal(lengths(ct$compounds$sample_idx_obs),
               lengths(ct$compounds$rt_frac_obs))
})

test_that("synthetic peaks land where padded real peaks are", {
  H <- 224L; W <- 48L
  smp <- make_toy_samples(8L, toy_peaks)
  ct <- build_peak_catalog(pp_from_samples(smp), eps_rt = 3, eps_cv = 1,
                           coord_mode = "physical")
  cm <- ct$compounds
  expect_equal(nrow(cm), nrow(toy_peaks))
  truth <- toy_peaks[order(toy_peaks$row), ]
  cm <- cm[order(cm$rt_loc), ]
  # where the real peaks sit after pad_to_target (12-row, 6-col offset)
  P <- pad_to_target(smp$Z[[1]], H, W)
  for (k in seq_len(nrow(truth)))
    expect_equal(P[12 + truth$row[k], 6 + truth$col[k]],
                 max(P[12 + truth$row[k] + -2:2, 6 + truth$col[k] + -1:1]))

  pos <- catalog_pixel_positions(ct, H, W)
  pos <- pos[match(cm$compound_id, pos$compound_id), ]
  expect_equal(pos$rt_px, 12 + truth$row)
  expect_equal(pos$cv_px, 6 + truth$col)

  res <- generate_one_synthetic_from_catalog(
    ct, H, W, add_noise = FALSE, size_jitter = 0,
    location_jitter_rt = 0, location_jitter_cv = 0,
    prevalence_override = stats::setNames(rep(1, nrow(cm)),
                                          as.character(cm$compound_id)))
  pl <- res$placements[match(cm$compound_id, res$placements$compound_id), ]
  expect_equal(pl$mu_rt, 12 + truth$row)
  expect_equal(pl$mu_cv, 6 + truth$col)
  # the old fraction * W rule put the high-CV peak ~5 columns too far right
  expect_gt(abs((truth$col[4] / 36) * W - (6 + truth$col[4])), 5)

  geo <- compound_head_geometry(ct, H, W)
  geo_px <- geo$cv_frac[match(cm$compound_id, geo$compound_id)] * W + 0.5
  expect_equal(geo_px, 6 + truth$col)
})

test_that("observations use their own sample's padding offset", {
  H <- 224L; W <- 48L
  a <- pp_from_samples(make_toy_samples(5L, toy_peaks, C = 36L, seed = 1L))
  b <- pp_from_samples(make_toy_samples(5L, toy_peaks, C = 40L, seed = 2L))
  ct <- build_universal_catalog(list(a = a, b = b), eps_rt = 3, eps_cv = 1,
                                min_support_frac = 0.1,
                                coord_mode = "physical")
  cm <- ct$compounds
  pl <- .catalog_placement(ct, cm, H, W, "test")
  i <- which.max(cm$cv_loc)                 # true col 34 in both cohorts
  for (j in seq_along(cm$rt_frac_obs[[i]])) {
    s <- cm$sample_idx_obs[[i]][j]
    off <- if (s <= 5L) 6 else 4            # floor((48 - C) / 2)
    expect_equal(.place_obs(pl, cm, i, j)[2], off + 34)
  }
})

test_that("dims are inferred from axis ranges; legacy catalogs scale with a warning", {
  smp <- make_toy_samples(4L, toy_peaks)
  pp <- pp_from_samples(smp)
  pp_old <- pp; pp_old$data_dims <- NULL
  ct_new <- build_peak_catalog(pp, eps_rt = 3, eps_cv = 1, coord_mode = "physical")
  rm(list = ls(.classydms_state), envir = .classydms_state)
  expect_message(
    ct_inf <- build_peak_catalog(pp_old, eps_rt = 3, eps_cv = 1,
                                 coord_mode = "physical"),
    "axis lengths")
  expect_equal(unname(ct_inf$parameters$data_dims), unname(pp$data_dims))
  expect_equal(catalog_pixel_positions(ct_inf, 224, 48),
               catalog_pixel_positions(ct_new, 224, 48))

  legacy <- ct_new
  legacy$parameters$data_dims <- NULL
  legacy$compounds$sample_idx_obs <- NULL
  rm(list = ls(.classydms_state), envir = .classydms_state)
  expect_warning(pos <- catalog_pixel_positions(legacy, 224, 48), "data dims")
  expect_equal(pos$cv_px, vapply(legacy$compounds$cv_frac_obs,
                                 function(v) stats::median(v * 48), numeric(1)))
})

test_that("peak_params generators and holdout templates use padded placement", {
  pp <- pp_from_samples(make_toy_samples(4L, toy_peaks))
  px <- .place_pp_peaks(pp, seq_along(pp$rt_loc_raw), 224, 48)
  expect_true(all(px[, 2] %in% (6 + toy_peaks$col)))
  expect_true(all(px[, 1] %in% (12 + toy_peaks$row)))
  ho <- holdout_biomarkers(pp, n_biomarkers = 1L, seed = 1L)
  expect_equal(ho$biomarker_templates$data_cols, 36L)
  pl <- list(H = 224, W = 48, ref_rows = 200, ref_cols = 36)
  expect_equal(.data_band_range(pl, c(0, 1), "cv"), c(6, 42) / 48)
})

test_that("pretrain samples keep contamination inside the data block", {
  H <- 224L; W <- 48L
  ct <- build_peak_catalog(pp_from_samples(make_toy_samples(6L, toy_peaks)),
                           eps_rt = 3, eps_cv = 1, coord_mode = "physical")
  set.seed(3)
  res <- generate_noisy_pretrain_sample(
    ct, H, W, anchor_ids = ct$compounds$compound_id,
    contamination_config = list(n_per_sample = 20, exclusion_eps_rt = 0.01,
                                exclusion_eps_cv = 0.02),
    add_noise = FALSE, location_jitter_rt = 0, location_jitter_cv = 0)
  cp <- res$contamination_positions
  expect_gt(nrow(cp), 0)
  # data block occupies rows 13..212, cols 7..42 on this canvas
  expect_true(all(cp$mu_cv >= 7 & cp$mu_cv <= 42))
  expect_true(all(cp$mu_rt >= 13 & cp$mu_rt <= 212))
})
