# =============================================================================
# Real-scale sanity-max on a peak catalog — end-to-end pipeline validation
# =============================================================================
#
# Runs:
#   1. Load an existing peak_params object (from May .RData or fresh)
#   2. (optional) backfill sample_idx_raw for May-era objects
#   3. (optional) subset by class for per-condition catalog
#   4. Build the peak catalog
#   5. Report catalog diagnostics
#   6. Simulate a case/control cohort using sanity-max biomarkers
#   7. Pretrain the encoder on catalog-based synthetic data
#   8. Encoder-channel Wilcoxon (per-channel class differentiation check)
#   9. Elastic net classification (5-fold CV)
#  10. Aggregate saliency + biomarker recovery evaluation
#  11. Save intermediate objects so you can rerun downstream steps
#      without redoing the expensive ones.
#
# Every step saves an .rds so a crash mid-run doesn't lose everything.
# The RUN_STEPS block at the top lets you skip steps you've already
# completed by setting them to FALSE (the previous save is loaded
# instead).
#
# Copy this file to your working directory and edit the CONFIG block
# below for your paths and hyperparameters.
# =============================================================================

library(classydms)
library(dplyr)
library(tibble)
library(purrr)

# =============================================================================
# CONFIG — edit these
# =============================================================================

# Path to your saved peak_params object (.RData or .rds)
PEAK_PARAMS_PATH <- "path/to/your/aki_peak_params.RData"

# If your saved file has a variable called "pp" (or whatever), tell us here.
# For .RData: assign after load(); for .rds: readRDS returns the value directly.
PEAK_PARAMS_VAR_NAME <- "pp"   # only used for .RData; ignored for .rds

# If your peak_params combines disease + controls, provide the class labels
# to subset by. Set to NULL to use all samples together.
# Class vector should be same length as n_samples in peak_params.
# Convention: 2 = disease (case), 1 = control.
CLASS_LABELS <- NULL             # e.g. c(rep(2, 88), rep(1, 101))
KEEP_CLASS  <- NULL              # e.g. c(1, 2) to keep both; c(1) for controls only

# Sample names (if you have them and want them backfilled).
# Should be a character vector of length n_samples in your peak_params.
SAMPLE_NAMES <- NULL

# Image dimensions (padded / trimmed sizes used for classification).
# Use whatever your real pipeline uses — for AKI native this is often
# 1400 x 40 or similar.
H <- 1400L
W <- 40L

# Catalog clustering parameters
EPS_RT <- 0.003
EPS_CV <- 0.015
MIN_SUPPORT_FRAC <- 0.10

# Case/control simulation
N_CASES    <- 100L
N_CONTROLS <- 100L

# Pretraining hyperparameters (production-ish, tuned down slightly for
# reasonable runtime on CPU — bump on GPU/MPS)
STEPS_PER_EPOCH   <- 300L
EPOCHS            <- 20L
BATCH_SIZE        <- 32L
LR                <- 1e-3
NOISE_SCALE       <- 3.0
SIZE_JITTER       <- 0.15
PRETRAIN_DEVICE   <- "cpu"   # "cuda" or "mps" if available
PRETRAIN_WORKERS  <- 1L      # bump to N-1 on Unix for parallel data gen

# Classification
CV_K     <- 5L
CV_ALPHA <- 0.5

# Where to save intermediate .rds objects
OUTPUT_DIR <- "catalog_sanity_max_run"

# Which steps to run. Set FALSE to skip a step and load the previous
# save instead. Useful for iterating on downstream steps after
# pretraining has completed.
RUN_STEPS <- list(
  load_peak_params = TRUE,
  build_catalog    = TRUE,
  simulate_cohort  = TRUE,
  pretrain_encoder = TRUE,
  extract_features = TRUE,
  classify         = TRUE,
  saliency         = TRUE
)

# Reproducibility
GLOBAL_SEED <- 42L

# =============================================================================
# END CONFIG — don't edit below unless you know what you're doing
# =============================================================================

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
saved_path <- function(name) file.path(OUTPUT_DIR, paste0(name, ".rds"))

save_step <- function(obj, name) {
  path <- saved_path(name)
  saveRDS(obj, path)
  cat(sprintf("  saved: %s\n", path))
}
load_step <- function(name) {
  path <- saved_path(name)
  if (!file.exists(path))
    stop("Missing saved step: ", path, " -- rerun that step (set to TRUE)")
  cat(sprintf("  loaded: %s\n", path))
  readRDS(path)
}

# =============================================================================
# STEP 1 — Load peak_params
# =============================================================================
cat("\n=== STEP 1: Load peak_params ===\n")
if (RUN_STEPS$load_peak_params) {
  if (grepl("\\.rds$", PEAK_PARAMS_PATH, ignore.case = TRUE)) {
    pp <- readRDS(PEAK_PARAMS_PATH)
  } else {
    # .RData: load into an environment, then extract the requested var
    env <- new.env()
    load(PEAK_PARAMS_PATH, envir = env)
    if (!PEAK_PARAMS_VAR_NAME %in% ls(env)) {
      stop(sprintf("Variable '%s' not found in %s. Available: %s",
                    PEAK_PARAMS_VAR_NAME, PEAK_PARAMS_PATH,
                    paste(ls(env), collapse = ", ")))
    }
    pp <- get(PEAK_PARAMS_VAR_NAME, envir = env)
  }
  cat(sprintf("  Loaded peak_params: %d samples, %d pool peaks\n",
              pp$n_samples, pp$n_peaks_detected))

  # Backfill sample_idx_raw if needed
  if (is.null(pp$sample_idx_raw)) {
    cat("  sample_idx_raw missing — backfilling...\n")
    pp <- backfill_sample_idx_raw(pp, sample_names = SAMPLE_NAMES)
    cat(sprintf("  backfilled sample_idx_raw (length %d)\n",
                length(pp$sample_idx_raw)))
  } else {
    cat("  sample_idx_raw already present\n")
  }

  # Subset by class if requested
  if (!is.null(CLASS_LABELS) && !is.null(KEEP_CLASS)) {
    keep_samples <- which(CLASS_LABELS %in% KEEP_CLASS)
    cat(sprintf("  Subsetting to %d samples in classes %s\n",
                length(keep_samples), paste(KEEP_CLASS, collapse = ",")))
    pp <- subset_peak_params(pp, keep_samples)
  }
  save_step(pp, "01_peak_params")
} else {
  pp <- load_step("01_peak_params")
}

# =============================================================================
# STEP 2 — Build catalog
# =============================================================================
cat("\n=== STEP 2: Build catalog ===\n")
if (RUN_STEPS$build_catalog) {
  t0 <- proc.time()[3]
  catalog <- build_peak_catalog(
    pp,
    eps_rt = EPS_RT, eps_cv = EPS_CV,
    min_support_frac = MIN_SUPPORT_FRAC
  )
  cat(sprintf("  Built in %.1f s\n", proc.time()[3] - t0))
  cat("\n")
  catalog_summary(catalog)
  cat("\nCompactness diagnostic:\n")
  print(catalog_compactness(catalog, H = H, W = W))
  save_step(catalog, "02_catalog")
} else {
  catalog <- load_step("02_catalog")
}

# =============================================================================
# STEP 3 — Simulate case/control cohort (sanity-max scenario)
# =============================================================================
cat("\n=== STEP 3: Simulate case/control cohort ===\n")
if (RUN_STEPS$simulate_cohort) {
  set.seed(GLOBAL_SEED)
  biom_spec <- catalog_biomarkers_sanity_max(
    catalog,
    min_base_prevalence = 0.3,
    seed = GLOBAL_SEED
  )
  cat("  Biomarker spec:\n")
  print(biom_spec)

  sim <- simulate_case_control_from_catalog(
    catalog, biomarkers = biom_spec,
    n_cases = N_CASES, n_controls = N_CONTROLS,
    H = H, W = W, seed = GLOBAL_SEED + 1L
  )
  cat(sprintf("  Cohort: %d samples (%d cases, %d controls)\n",
              length(sim$samples), sum(sim$y == 2L), sum(sim$y == 1L)))

  # Verify: cases fire all biomarkers, controls fire none
  case_fire_pct <- mean(rowSums(sim$per_sample_firing[sim$y == 2L, ,
                                                        drop = FALSE]) ==
                          nrow(biom_spec))
  ctrl_fire_pct <- mean(rowSums(sim$per_sample_firing[sim$y == 1L, ,
                                                        drop = FALSE]) > 0)
  cat(sprintf("  Cases firing ALL biomarkers: %.0f%% (should be 100%%)\n",
              100 * case_fire_pct))
  cat(sprintf("  Controls firing ANY biomarker: %.0f%% (should be 0%%)\n",
              100 * ctrl_fire_pct))

  save_step(list(sim = sim, biom_spec = biom_spec), "03_simulation")
} else {
  loaded <- load_step("03_simulation")
  sim <- loaded$sim; biom_spec <- loaded$biom_spec
}

# =============================================================================
# STEP 4 — Pretrain encoder on catalog-based synthetic data
# =============================================================================
cat("\n=== STEP 4: Pretrain encoder ===\n")
if (RUN_STEPS$pretrain_encoder) {
  t0 <- proc.time()[3]
  pretrain_result <- pretrain_denoising_online(
    peak_params      = NULL,
    catalog          = catalog,
    H = H, W = W,
    steps_per_epoch  = STEPS_PER_EPOCH,
    epochs           = EPOCHS,
    batch_size       = BATCH_SIZE,
    lr               = LR,
    size_jitter      = SIZE_JITTER,
    noise_scale      = NOISE_SCALE,
    val_n            = 64L,
    device           = PRETRAIN_DEVICE,
    num_workers      = PRETRAIN_WORKERS,
    seed             = GLOBAL_SEED + 2L
  )
  cat(sprintf("  Pretrain complete in %.1f min\n",
              (proc.time()[3] - t0) / 60))
  cat(sprintf("  Final train MSE: %.4f, val MSE: %.4f\n",
              tail(pretrain_result$loss_history, 1),
              tail(pretrain_result$val_loss_history, 1)))
  save_step(pretrain_result, "04_pretrain")
} else {
  pretrain_result <- load_step("04_pretrain")
}
encoder <- pretrain_result$encoder

# =============================================================================
# STEP 5 — Extract features + per-channel Wilcoxon
# =============================================================================
cat("\n=== STEP 5: Encoder channel differentiation ===\n")
if (RUN_STEPS$extract_features) {
  Z_list <- lapply(sim$samples, `[[`, "Z")
  feats <- extract_encoder_features(encoder, Z_list, device = "cpu")

  p_raw <- sapply(seq_len(64L), function(c)
    suppressWarnings(wilcox.test(feats[sim$y == 2L, c],
                                  feats[sim$y == 1L, c])$p.value))
  p_bh <- p.adjust(p_raw, method = "BH")
  n_sig_05  <- sum(p_bh < 0.05, na.rm = TRUE)
  n_sig_001 <- sum(p_bh < 0.001, na.rm = TRUE)
  cat(sprintf("  Encoder channels significant at BH q < 0.05:  %d / 64\n",
              n_sig_05))
  cat(sprintf("  Encoder channels significant at BH q < 0.001: %d / 64\n",
              n_sig_001))

  save_step(list(features = feats, p_bh = p_bh,
                 n_sig_05 = n_sig_05, n_sig_001 = n_sig_001),
             "05_features")
} else {
  loaded <- load_step("05_features")
  feats <- loaded$features; p_bh <- loaded$p_bh
}

# =============================================================================
# STEP 6 — Classify with elastic net
# =============================================================================
cat("\n=== STEP 6: Classification ===\n")
if (RUN_STEPS$classify) {
  Z_list <- lapply(sim$samples, `[[`, "Z")
  elnet <- cv_encoder_elastic_net(
    encoder, Z_list, sim$y,
    k = CV_K, alpha = CV_ALPHA, seed = GLOBAL_SEED + 3L
  )
  cat(sprintf("  Case/control AUC (5-fold CV): %.3f\n", elnet$auc))
  cat(sprintf("  Non-zero elastic-net coefs: %d / 64\n",
              sum(elnet$coefs != 0)))
  save_step(elnet, "06_elnet")
} else {
  elnet <- load_step("06_elnet")
}

# =============================================================================
# STEP 7 — Aggregate saliency + biomarker recovery
# =============================================================================
cat("\n=== STEP 7: Saliency and biomarker recovery ===\n")
if (RUN_STEPS$saliency) {
  Z_list <- lapply(sim$samples, `[[`, "Z")
  case_indices <- which(sim$y == 2L)
  orig_dims <- lapply(sim$samples,
                       function(s) c(H = nrow(s$Z), W = ncol(s$Z)))

  sal <- aggregate_saliency_masked(encoder, Z_list,
                                     coefs = elnet$coefs,
                                     sample_indices = case_indices,
                                     orig_dims = orig_dims)
  # For biomarker recovery, we need the biomarker locations in the same
  # (H, W) coordinate frame as the saliency map.
  # biom_spec has compound_id; look them up in catalog for rt/cv.
  biom_locs <- catalog$compounds[
    match(biom_spec$compound_id, catalog$compounds$compound_id),
    c("rt_loc", "cv_loc")
  ]

  recovery <- evaluate_biomarker_recovery(
    sal, biom_locs,
    H = H, W = W,
    match_radius_rt = 8L, match_radius_cv = 2L,
    top_k = nrow(biom_spec)
  )
  cat(sprintf("  Biomarker recovery per-pixel AUC: %.3f\n",
              recovery$auc_per_pixel))
  cat(sprintf("  Recall @ top-%d: %.3f\n",
              nrow(biom_spec), recovery$recall_at_k))
  cat(sprintf("  Precision @ top-%d: %.3f\n",
              nrow(biom_spec), recovery$precision_at_k))
  save_step(list(sal = sal, recovery = recovery), "07_saliency")
} else {
  loaded <- load_step("07_saliency")
  sal <- loaded$sal; recovery <- loaded$recovery
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================
cat("\n\n=== FINAL SUMMARY ===\n")
cat(sprintf("Catalog:           %d compounds at eps=(%.3f, %.3f)\n",
            nrow(catalog$compounds), EPS_RT, EPS_CV))
cat(sprintf("Cohort:            %d cases, %d controls\n",
            N_CASES, N_CONTROLS))
cat(sprintf("Encoder channels:  %d / 64 significant at BH q<0.05\n",
            sum(p_bh < 0.05, na.rm = TRUE)))
cat(sprintf("Case/control AUC:  %.3f\n", elnet$auc))
cat(sprintf("Biomarker recovery AUC:  %.3f\n", recovery$auc_per_pixel))
cat(sprintf("Recall @ top-%d:   %.3f\n",
            nrow(biom_spec), recovery$recall_at_k))
cat("\n")
cat("Acceptance targets (sanity-max, catalog-based):\n")
cat("  Case/control AUC:      >= 0.95 (target reached: ",
    if (elnet$auc >= 0.95) "YES" else "NO", ")\n", sep = "")
cat("  Biomarker recovery:    >= 0.85 (target reached: ",
    if (recovery$auc_per_pixel >= 0.85) "YES" else "NO", ")\n", sep = "")
cat("  Encoder channels sig:  >= 30 / 64 (target reached: ",
    if (sum(p_bh < 0.05, na.rm = TRUE) >= 30L) "YES" else "NO", ")\n", sep = "")
cat("\n")
cat("All intermediates saved in ", OUTPUT_DIR, "\n", sep = "")
cat("To rerun a step: set RUN_STEPS$<step_name> = FALSE for steps you\n")
cat("want to skip (they'll be loaded from disk), TRUE for steps to rerun.\n")
