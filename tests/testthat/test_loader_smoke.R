# Smoke test against local sample data (not run during R CMD check —
# sampledata/ is excluded from the package build via .Rbuildignore).
# Run interactively from the package root to verify the full loader
# + preprocess + estimate_peak_params pipeline against real files.

if (interactive() &&
    dir.exists("sampledata/influenza") &&
    dir.exists("sampledata/controls")) {

  # Source the R files (so we don't have to install first)
  for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
    source(f)
  }

  loaded <- load_dms_directory(
    disease_dir = "sampledata/influenza",
    control_dir = "sampledata/controls",
    pattern = "POS\\.xls$"
  )

  stopifnot(
    is.list(loaded$samples),
    length(loaded$samples) > 0,
    is.null(names(loaded$samples)),  # no zero-length name risk
    is.matrix(loaded$samples[[1]]$Z),
    is.null(dimnames(loaded$samples[[1]]$Z))
  )

  samples_proc <- lapply(loaded$samples, process_one_sample,
                          sg_p = 3, sg_n = 21,
                          asls_lambda = 1e6, asls_p = 0.01, asls_niter = 10)

  samples_dust <- lapply(samples_proc, baseline_basement, basement_thr = 0.005)

  b_list <- lapply(samples_dust,
                    function(s) trim_bounds_from_occupancy(s$Z))
  R0 <- max(sapply(samples_dust, function(s) nrow(s$Z)))
  C0 <- max(sapply(samples_dust, function(s) ncol(s$Z)))
  global_b <- pooled_trim_bounds(b_list, R = R0, C = C0)
  samples_trimmed <- lapply(samples_dust, trim_sample, bounds = global_b)

  Z_norm <- lapply(samples_trimmed, function(s) normalize_sample(s$Z))
  pad_dims <- compute_pad_targets(Z_norm)
  Z_padded <- lapply(Z_norm,
                      function(z) pad_to_target(z, pad_dims$H, pad_dims$W))

  Z_raw <- lapply(samples_trimmed, function(s) s$Z)
  Z_pretrim <- lapply(samples_trimmed,
                       function(s) if (!is.null(s$Z_pretrim)) s$Z_pretrim else s$Z)

  pp <- estimate_peak_params(Z_raw, Z_pretrim_list = Z_pretrim,
                              top_k = 50, frac_height = 0.5)

  message("\nSmoke test passed:")
  message("  Loaded ", length(loaded$samples), " samples")
  message("  Pad targets: ", pad_dims$H, " x ", pad_dims$W)
  message("  Peak params detected: ", pp$n_peaks_detected, " peaks")
  message("  sigma_rt meanlog: ", round(pp$sigma_rt$meanlog, 2))
  message("  sigma_cv meanlog: ", round(pp$sigma_cv$meanlog, 2))
}
