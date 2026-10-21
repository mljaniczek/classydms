test_that("subset_peak_params keeps every per-peak vector aligned", {
  pp <- list(rt_loc_raw = c(.1, .2, .3, .4), cv_loc_raw = rep(.5, 4),
             sigma_rt_raw = rep(2, 4), sigma_cv_raw = rep(1, 4),
             intensity_raw = c(1, 2, 3, 4), sample_idx_raw = c(1L, 1L, 2L, 2L),
             rt_seconds_raw = c(150, 300, 450, 600), cv_volts_raw = rep(-5, 4),
             n_samples = 2L, sample_names = c("a", "b"),
             coord_mode = "physical")
  s <- subset_peak_params(pp, 2L)
  expect_equal(s$rt_loc_raw, c(.3, .4))
  expect_equal(s$rt_seconds_raw, c(450, 600))
  expect_equal(s$cv_volts_raw, c(-5, -5))
  expect_equal(s$sample_idx_raw, c(1L, 1L))
  expect_equal(s$sample_names, "b")
  expect_silent(validate_peak_params(s))
})

test_that("subset_peak_params subsets per-sample fields", {
  pp <- make_toy_peak_params(6L, toy_compounds())
  pp$time_ranges[, "max"] <- 100 + seq_len(6)
  pp$alignment <- list(shifts_rt = seq_len(6) / 10, shifts_cv = -seq_len(6),
                       n_anchors_matched = 10L + seq_len(6),
                       n_anchors_total = 4L)
  s <- subset_peak_params(pp, c(5L, 2L))
  expect_equal(s$n_samples, 2L)
  expect_equal(unname(s$time_ranges[, "max"]), c(102, 105))
  expect_equal(nrow(s$cv_ranges), 2L)
  expect_equal(s$n_peaks$values, pp$n_peaks$values[c(2, 5)])
  expect_equal(s$alignment$shifts_rt, c(.2, .5))
  expect_equal(s$alignment$n_anchors_total, 4L)
  expect_equal(s$sample_names, pp$sample_names[c(2, 5)])
  expect_error(subset_peak_params(pp, c(2L, 2L)), "duplicate")
})

test_that("validate_peak_params rejects misaligned objects", {
  pp <- make_toy_peak_params(4L, toy_compounds())
  expect_silent(validate_peak_params(pp))
  bad <- pp; bad$rt_seconds_raw <- bad$rt_seconds_raw[-1]
  expect_error(validate_peak_params(bad), "rt_seconds_raw")
  bad <- pp; bad$sample_idx_raw[1] <- 99L
  expect_error(validate_peak_params(bad), "sample_idx_raw")
  bad <- pp; bad$time_ranges <- bad$time_ranges[-1, ]
  expect_error(validate_peak_params(bad), "time_ranges")
  bad <- pp; bad$cv_volts_raw <- NULL
  expect_error(validate_peak_params(bad), "cv_volts_raw")
  # build functions refuse a misaligned object instead of NA-dropping
  bad <- pp; bad$rt_loc_raw <- bad$rt_loc_raw[-1]
  expect_error(build_peak_catalog(bad, eps_rt = 3, eps_cv = 1,
                                  coord_mode = "physical"), "rt_loc_raw")
})

test_that("merge_peak_params carries physical and per-sample fields", {
  a <- make_toy_peak_params(5L, toy_compounds(), seed = 1L)
  b <- make_toy_peak_params(7L, toy_compounds(), seed = 3L)
  m <- merge_peak_params(list(a = a, b = b))
  expect_silent(validate_peak_params(m))
  expect_equal(m$rt_seconds_raw, c(a$rt_seconds_raw, b$rt_seconds_raw))
  expect_equal(m$cv_volts_raw, c(a$cv_volts_raw, b$cv_volts_raw))
  expect_equal(nrow(m$time_ranges), 12L)
  expect_equal(m$sample_names, c(a$sample_names, b$sample_names))
  expect_equal(m$coord_mode, "physical")
  # physical fields dropped (with a message) when not all inputs have them
  b2 <- b; b2$rt_seconds_raw <- NULL; b2$cv_volts_raw <- NULL
  expect_message(m2 <- merge_peak_params(list(a = a, b = b2)), "rt_seconds_raw")
  expect_null(m2$rt_seconds_raw)
  expect_silent(validate_peak_params(m2))
})

# End-to-end: subset (drop holdout samples) -> merge -> physical build.
# This is the path that produced scrambled catalogs before the fix.
test_that("physical universal catalog from subsetted cohorts is consistent", {
  R <- 300L; C <- 40L
  cmp <- toy_compounds(R = R, C = C)
  # Samples 1:4 of cohort a are the holdout. They alone carry a compound
  # at row 150 / col 3 that must not appear in the catalog.
  priv <- list(samples = 1:4,
               compounds = data.frame(row = 150, col = 3, intensity = 50,
                                      sigma_rt = 99, prevalence = 1))
  a <- make_toy_peak_params(30L, cmp, R = R, C = C, private = priv, seed = 11L)
  b <- make_toy_peak_params(25L, cmp, R = R, C = C, seed = 12L)
  a_kept <- subset_peak_params(a, 5:30)
  pp_list <- list(a = a_kept, b = b)
  ct <- build_universal_catalog(pp_list, eps_rt = 3, eps_cv = 1,
                                min_support_frac = 0.05,
                                coord_mode = "physical")
  merged <- merge_peak_params(pp_list)
  cm <- ct$compounds
  expect_equal(nrow(cm), nrow(cmp))
  expect_equal(ct$parameters$n_pool_peaks, length(merged$rt_loc_raw))
  expect_equal(ct$parameters$n_samples, 51L)

  # 1. Positions: every observation's fraction maps back to its compound.
  for (i in seq_len(nrow(cm))) {
    rows_obs <- cm$rt_frac_obs[[i]] * R
    secs_row <- cm$rt_loc[i] + 1          # time = row - 1 on this grid
    expect_true(all(abs(rows_obs - secs_row) <= 2))
    cols_obs <- cm$cv_frac_obs[[i]] * C
    expect_true(all(abs(cols_obs - ((cm$cv_loc[i] + 20) / 0.5 + 1)) <= 0.5))
  }
  expect_gt(stats::cor(cm$rt_frac, cm$rt_loc), 0.999)
  expect_gt(stats::cor(cm$cv_frac, cm$cv_loc), 0.999)

  # 2. Morphology/intensity belong to the compound's own peaks.
  truth <- vapply(cm$rt_loc, function(s) which.min(abs(cmp$row - (s + 1))),
                  integer(1))
  for (i in seq_len(nrow(cm))) {
    expect_true(all(cm$sigma_rt_obs[[i]] == cmp$sigma_rt[truth[i]]))
    expect_true(all(abs(log(cm$intensity_obs[[i]] / cmp$intensity[truth[i]])) < 0.3))
    expect_equal(lengths(cm$rt_frac_obs)[i], cm$n_observations[i])
  }

  # 3. Excluded samples contribute nothing.
  expect_false(any(abs(cm$rt_loc - 149) < 5 & abs(cm$cv_loc - (-19)) < 2))
  expect_false(any(unlist(cm$sigma_rt_obs) == 99))
  all_members <- unlist(cm$member_indices)
  expect_true(all(all_members <= length(merged$rt_loc_raw)))

  # 4. Prevalence = unique contributing samples / n_samples, per source too.
  for (i in seq_len(nrow(cm))) {
    s <- unique(merged$sample_idx_raw[cm$member_indices[[i]]])
    expect_equal(cm$n_samples[i], length(s))
    expect_equal(cm$prevalence[i], length(s) / 51)
    src <- table(factor(merged$sample_cohort[s], levels = c("a", "b")))
    expect_equal(unname(cm$per_source_prevalence[[i]]),
                 as.numeric(src) / c(26, 25))
  }
  expect_true(all(cm$prevalence <= 1))
})

test_that("per-arm physical catalog from a subset is consistent", {
  R <- 300L; C <- 40L
  a <- make_toy_peak_params(40L, toy_compounds(R = R, C = C), R = R, C = C,
                            seed = 5L)
  arm <- subset_peak_params(a, seq(2L, 40L, by = 2L))
  ct <- build_peak_catalog(arm, eps_rt = 3, eps_cv = 1,
                           min_support_frac = 0.05, coord_mode = "physical")
  cm <- ct$compounds
  expect_gt(stats::cor(cm$rt_frac * R, cm$rt_loc), 0.999)
  for (i in seq_len(nrow(cm))) {
    expect_true(all(abs(cm$rt_frac_obs[[i]] * R - (cm$rt_loc[i] + 1)) <= 2))
    s <- unique(arm$sample_idx_raw[cm$member_indices[[i]]])
    expect_equal(cm$prevalence[i], length(s) / 20)
  }
  # fraction-mode catalog's rt_seconds also come from the right peaks
  ctf <- build_peak_catalog(arm, eps_rt = 3 / R, eps_cv = 1 / C,
                            min_support_frac = 0.05, coord_mode = "fraction")
  expect_gt(stats::cor(ctf$compounds$rt_frac * R,
                       ctf$compounds$rt_seconds), 0.999)
})

test_that("holdout_biomarkers keeps background vectors aligned", {
  a <- make_toy_peak_params(10L, toy_compounds())
  h <- holdout_biomarkers(a, n_biomarkers = 2L, seed = 1L)
  expect_silent(validate_peak_params(h$background_peak_params))
})

test_that("alignment paths keep peak_params aligned and sync physical", {
  cmp <- toy_compounds()
  a <- make_toy_peak_params(20L, cmp, seed = 7L)
  b <- make_toy_peak_params(20L, cmp, seed = 8L)
  # shift cohort b by 4 rows (4 s) so alignment has something to do
  b$rt_loc_raw <- b$rt_loc_raw + 4 / 300
  b$rt_seconds_raw <- b$rt_seconds_raw + 4
  al <- align_across_cohorts(list(a = subset_peak_params(a, 3:20), b = b),
                             eps_rt_loose = 0.05, eps_cv_loose = 0.1,
                             min_samples_per_anchor = 5L, verbose = FALSE)
  for (pp in al) expect_silent(validate_peak_params(pp))
  # physical coordinates moved with the fractions
  d_frac <- (al$b$rt_loc_raw - b$rt_loc_raw) * 300
  d_sec  <- al$b$rt_seconds_raw - b$rt_seconds_raw
  expect_lt(max(abs(d_frac - d_sec)), 0.1)
  expect_gt(abs(stats::median(d_sec)), 1)

  ap <- suppressWarnings(align_peak_params(a, verbose = FALSE))
  expect_silent(validate_peak_params(ap))
  expect_true(ap$alignment$physical_synced)

  ct <- suppressMessages(build_universal_catalog(
    list(a = a, b = b), eps_rt = 3, eps_cv = 1, min_support_frac = 0.05,
    coord_mode = "physical", align = list(min_samples_per_anchor = 5L)))
  expect_gt(stats::cor(ct$compounds$rt_frac * 300, ct$compounds$rt_loc), 0.999)
})
