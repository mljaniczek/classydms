# classydms: case / control simulation harness
# Ground-truth validation of the full pipeline. Simulates a synthetic
# cohort where cases differ from controls by known injected peaks
# ("biomarkers") with configurable prevalence and intensity differences,
# then feeds through preprocess -> params -> pretrain -> classify ->
# saliency to check that the pipeline recovers the injected signal.

#' Hold out K observed peaks as biomarker templates
#'
#' Uniformly at random draws `n_biomarkers` peak indices from
#' `peak_params`'s raw pools, extracts their locations and morphology
#' as biomarker templates, and returns a modified `peak_params` with
#' those peaks *and* any peaks within a protected radius removed from
#' the raw pools. This guarantees that background samples drawn from
#' the returned `background_peak_params` will not include peaks that
#' land near biomarker locations, so injected biomarkers are truly
#' additive signal rather than contaminated by nearby background peaks.
#'
#' @param peak_params Output of [estimate_peak_params()]; must contain
#'   `rt_loc_raw`, `cv_loc_raw`, `sigma_rt_raw`, `sigma_cv_raw`,
#'   `intensity_raw`.
#' @param n_biomarkers Integer; number of biomarker templates to draw.
#' @param protect_radius_rt,protect_radius_cv Radii (in fractions of
#'   image dimensions) around each biomarker inside which observed
#'   peaks are also removed from the background pool. Default matches
#'   the default `min_sep_rt`, `min_sep_cv` on a 1400 x 40 image (15 /
#'   1400 and 3 / 40).
#' @param min_intensity_quantile Optional numeric in `[0, 1)`. When
#'   non-`NULL`, biomarker templates are drawn only from pool peaks
#'   whose `intensity_raw` is at or above the given quantile of the
#'   pool's intensity distribution — e.g. `0.90` restricts selection
#'   to the top-10% brightest peaks. Useful for sanity-check
#'   scenarios where you want the injected biomarker to be visually
#'   obvious even in a very sparse image. `NULL` (default) samples
#'   uniformly from the whole pool.
#' @param seed Optional RNG seed for reproducible template selection.
#'
#' @return A list with:
#' \describe{
#'   \item{biomarker_templates}{tibble with columns `rt_loc`, `cv_loc`,
#'     `sigma_rt`, `sigma_cv`, `intensity` — one row per held-out peak.}
#'   \item{background_peak_params}{modified `peak_params` with the
#'     held-out peaks and any peaks within the protected radius stripped
#'     from every `*_raw` vector. `n_peaks_detected` is updated.}
#'   \item{n_pool_before}{pool size before hold-out.}
#'   \item{n_pool_after}{pool size after hold-out.}
#' }
#' @export
holdout_biomarkers <- function(peak_params,
                                n_biomarkers,
                                protect_radius_rt = 15 / 1400,
                                protect_radius_cv = 3  / 40,
                                min_intensity_quantile = NULL,
                                seed = NULL) {
  needed <- c("rt_loc_raw", "cv_loc_raw",
              "sigma_rt_raw", "sigma_cv_raw", "intensity_raw")
  missing <- needed[vapply(needed,
                            function(nm) is.null(peak_params[[nm]]),
                            logical(1))]
  if (length(missing)) {
    stop("holdout_biomarkers: peak_params is missing ",
         paste(missing, collapse = ", "), ".")
  }
  n_pool <- length(peak_params$rt_loc_raw)
  if (n_biomarkers >= n_pool) {
    stop("n_biomarkers (", n_biomarkers,
         ") must be less than pool size (", n_pool, ").")
  }
  if (!is.null(min_intensity_quantile)) {
    if (min_intensity_quantile < 0 || min_intensity_quantile >= 1) {
      stop("min_intensity_quantile must be in [0, 1); got ",
           min_intensity_quantile)
    }
  }
  if (!is.null(seed)) set.seed(seed)

  # Restrict the eligible pool by intensity, if requested
  eligible <- seq_len(n_pool)
  if (!is.null(min_intensity_quantile)) {
    thr <- as.numeric(stats::quantile(peak_params$intensity_raw,
                                        min_intensity_quantile))
    eligible <- which(peak_params$intensity_raw >= thr)
    if (length(eligible) < n_biomarkers) {
      stop("Only ", length(eligible), " peaks at or above the ",
           min_intensity_quantile, " intensity quantile, but ",
           n_biomarkers, " biomarkers requested.")
    }
  }
  biomarker_idx <- if (length(eligible) == n_biomarkers) eligible
                   else sample(eligible, size = n_biomarkers, replace = FALSE)

  templates <- tibble::tibble(
    rt_loc    = peak_params$rt_loc_raw[biomarker_idx],
    cv_loc    = peak_params$cv_loc_raw[biomarker_idx],
    sigma_rt  = peak_params$sigma_rt_raw[biomarker_idx],
    sigma_cv  = peak_params$sigma_cv_raw[biomarker_idx],
    intensity = peak_params$intensity_raw[biomarker_idx]
  )

  # Find every pool peak within the protected radius of ANY biomarker
  protected_mask <- logical(n_pool)   # FALSE by default
  for (b in seq_len(n_biomarkers)) {
    dr <- abs(peak_params$rt_loc_raw - templates$rt_loc[b])
    dc <- abs(peak_params$cv_loc_raw - templates$cv_loc[b])
    protected_mask <- protected_mask |
      (dr <= protect_radius_rt & dc <= protect_radius_cv)
  }
  # Also drop the biomarker peaks themselves (some are already inside
  # their own radius, but a zero-radius edge case would miss them)
  protected_mask[biomarker_idx] <- TRUE
  keep_mask <- !protected_mask

  bg <- peak_params
  for (nm in c("rt_loc_raw", "cv_loc_raw", "sigma_rt_raw",
               "sigma_cv_raw", "intensity_raw")) {
    bg[[nm]] <- peak_params[[nm]][keep_mask]
  }
  if (!is.null(bg$sample_idx_raw))
    bg$sample_idx_raw <- peak_params$sample_idx_raw[keep_mask]
  bg$n_peaks_detected <- sum(keep_mask)

  list(
    biomarker_templates    = templates,
    background_peak_params = bg,
    n_pool_before          = n_pool,
    n_pool_after           = sum(keep_mask)
  )
}

#' Combine biomarker templates with case/control differential spec
#'
#' Wraps templates in the group-level behavior spec used by
#' [simulate_case_control_cohort()]. Prevalences and intensity
#' multipliers may be scalars (broadcast to every biomarker) or
#' vectors of length matching `nrow(templates)`.
#'
#' @param templates tibble returned by [holdout_biomarkers()].
#' @param case_prevalence,control_prevalence Bernoulli probability that
#'   each biomarker appears in a case or control sample. Scalar or
#'   length-matched vector.
#' @param case_intensity_mult,control_intensity_mult Multiplier applied
#'   to the template `intensity` when the biomarker fires in a case or
#'   control sample. Scalar or length-matched vector.
#'
#' @return tibble with the template columns plus
#'   `case_prevalence`, `control_prevalence`, `case_intensity_mult`,
#'   `control_intensity_mult`. Suitable for `biomarkers = ...` in
#'   [simulate_case_control_cohort()].
#' @export
build_biomarker_spec <- function(templates,
                                  case_prevalence,
                                  control_prevalence,
                                  case_intensity_mult    = 1.0,
                                  control_intensity_mult = 1.0) {
  n <- nrow(templates)
  broadcast <- function(x, nm) {
    if (length(x) == 1L) return(rep(x, n))
    if (length(x) == n)  return(x)
    stop(nm, " must have length 1 or ", n, ", got ", length(x), ".")
  }
  templates$case_prevalence        <- broadcast(case_prevalence,
                                                 "case_prevalence")
  templates$control_prevalence     <- broadcast(control_prevalence,
                                                 "control_prevalence")
  templates$case_intensity_mult    <- broadcast(case_intensity_mult,
                                                 "case_intensity_mult")
  templates$control_intensity_mult <- broadcast(control_intensity_mult,
                                                 "control_intensity_mult")
  templates
}

#' Maximally-favorable ("sanity-max") biomarker scenario
#'
#' Combines every lever we have in the direction of "make signal
#' impossible to miss":
#' \itemize{
#'   \item 5 biomarkers (collective signal, not one lone peak)
#'   \item Drawn from top 5% of pool by intensity
#'   \item Perfect presence / absence: `case_prevalence = 1.0`,
#'     `control_prevalence = 0.0`
#'   \item Case intensity multiplier = 10x
#'   \item Sparse background: `n_peaks$mean = 30` and
#'     `n_peaks$sd = 5` overwritten on the returned peak_params so
#'     synthetic samples have only ~30 background peaks rather than
#'     ~1300, making injected biomarkers dominate total image intensity
#' }
#'
#' Intended as a **pipeline validation** scenario, not a difficulty
#' benchmark: if case/control AUC on this scenario is < 0.95, something
#' in the pipeline is broken and Levels 1/2/3 are not worth running.
#'
#' @inheritParams biomarkers_level1
#' @return list(biomarkers, background_peak_params) — pass both to
#'   [simulate_case_control_cohort()].
#' @export
biomarkers_sanity_max <- function(peak_params, seed = NULL) {
  ho <- holdout_biomarkers(peak_params, n_biomarkers = 5L,
                            min_intensity_quantile = 0.95,
                            seed = seed)
  bg <- ho$background_peak_params
  # Force sparse background so injected biomarkers dominate.
  bg$n_peaks$mean <- 30
  bg$n_peaks$sd   <- 5
  spec <- build_biomarker_spec(
    ho$biomarker_templates,
    case_prevalence        = 1.0,
    control_prevalence     = 0.0,
    case_intensity_mult    = 10.0,
    control_intensity_mult = 1.0
  )
  list(biomarkers = spec, background_peak_params = bg)
}

#' Level 1 (sanity check) biomarker scenario
#'
#' 1 biomarker drawn from the *top 10% brightest* pool peaks
#' (`min_intensity_quantile = 0.90`), extreme prevalence gap (95%
#' cases vs 5% controls), 5x intensity multiplier in cases. Designed
#' to be visually obvious in a spot-check plot; case/control AUC
#' should be well above 0.90 if the pipeline works. If Level 1 fails,
#' something in the pipeline is broken.
#'
#' @param peak_params Output of [estimate_peak_params()].
#' @param seed Optional RNG seed for reproducible biomarker selection.
#' @return list(biomarkers, background_peak_params) — pass both to
#'   [simulate_case_control_cohort()].
#' @export
biomarkers_level1 <- function(peak_params, seed = NULL) {
  ho <- holdout_biomarkers(peak_params, n_biomarkers = 1L,
                            min_intensity_quantile = 0.90,
                            seed = seed)
  spec <- build_biomarker_spec(
    ho$biomarker_templates,
    case_prevalence        = 0.95,
    control_prevalence     = 0.05,
    case_intensity_mult    = 5.0,
    control_intensity_mult = 1.0
  )
  list(biomarkers = spec, background_peak_params = ho$background_peak_params)
}

#' Level 2 (realistic) biomarker scenario
#'
#' 3 biomarkers drawn from the *upper half* of the pool by intensity
#' (`min_intensity_quantile = 0.50`), moderate prevalence gap
#' (~80% vs ~20%), 1.5x intensity multiplier in cases. Case/control
#' AUC should be at least ~0.85 if the pipeline works.
#' @inheritParams biomarkers_level1
#' @export
biomarkers_level2 <- function(peak_params, seed = NULL) {
  ho <- holdout_biomarkers(peak_params, n_biomarkers = 3L,
                            min_intensity_quantile = 0.50,
                            seed = seed)
  spec <- build_biomarker_spec(
    ho$biomarker_templates,
    case_prevalence        = 0.80,
    control_prevalence     = 0.20,
    case_intensity_mult    = 1.5,
    control_intensity_mult = 1.0
  )
  list(biomarkers = spec, background_peak_params = ho$background_peak_params)
}

#' Level 3 (hard) biomarker scenario
#'
#' 6 biomarkers drawn uniformly from the pool (any intensity), small
#' prevalence gap (~60% vs ~40%), 1.2x intensity multiplier in cases.
#' Stress-tests the pipeline including dim biomarkers; AUC around
#' 0.70 is a passing result at this difficulty.
#' @inheritParams biomarkers_level1
#' @export
biomarkers_level3 <- function(peak_params, seed = NULL) {
  ho <- holdout_biomarkers(peak_params, n_biomarkers = 6L,
                            min_intensity_quantile = NULL,
                            seed = seed)
  spec <- build_biomarker_spec(
    ho$biomarker_templates,
    case_prevalence        = 0.60,
    control_prevalence     = 0.40,
    case_intensity_mult    = 1.2,
    control_intensity_mult = 1.0
  )
  list(biomarkers = spec, background_peak_params = ho$background_peak_params)
}

#' Simulate a case / control cohort with known injected biomarkers
#'
#' Generates a synthetic cohort of `n_cases` case samples and
#' `n_controls` control samples. Each sample gets:
#' \enumerate{
#'   \item A background drawn from `background_peak_params` via
#'     [generate_one_synthetic()] with default `location_mode =
#'     "empirical"` and `attribute_mode = "joint"` (add_noise = FALSE
#'     at this stage; noise is added after biomarkers).
#'   \item For each biomarker row: a Bernoulli(`case_prevalence` or
#'     `control_prevalence`) draw decides whether the biomarker
#'     "fires" in this sample. If it fires, a 2D Gaussian peak with
#'     the biomarker's `(rt_loc, cv_loc, sigma_rt, sigma_cv)` and
#'     intensity `intensity * (case_intensity_mult | control_intensity_mult)`
#'     is added, plus per-peak jitter matching production
#'     (`location_jitter_rt/cv`, `size_jitter`).
#'   \item Folded-Normal background noise (from
#'     `background_peak_params$noise`).
#' }
#'
#' @param background_peak_params Modified `peak_params` returned by
#'   [holdout_biomarkers()] (or the same, minus the biomarker peaks).
#' @param biomarkers Output of [build_biomarker_spec()] (or a compatible
#'   tibble with the same columns).
#' @param n_cases,n_controls Number of samples in each group.
#' @param H,W Output image dimensions.
#' @param size_jitter,location_jitter_rt,location_jitter_cv Passed to
#'   the background generator; also applied to biomarker peaks so their
#'   morphology varies per sample.
#' @param add_noise Whether to add background noise (default TRUE).
#' @param seed Optional RNG seed.
#'
#' @return A list with:
#' \describe{
#'   \item{samples}{named list of length `n_cases + n_controls`. Each
#'     element is `list(path, time, cv, Z)` matching the shape produced
#'     by [load_dms_directory()], so it can flow through the rest of
#'     the classydms pipeline unchanged.}
#'   \item{y}{named integer vector: 2 = case, 1 = control.}
#'   \item{ground_truth}{the biomarkers tibble.}
#'   \item{per_sample_firing}{logical matrix (`n_samples` x
#'     `n_biomarkers`) recording which biomarkers actually fired in
#'     each sample, for post-hoc analysis.}
#' }
#' @export
simulate_case_control_cohort <- function(background_peak_params,
                                          biomarkers,
                                          n_cases, n_controls,
                                          H, W,
                                          size_jitter        = 0.15,
                                          location_jitter_rt = 2,
                                          location_jitter_cv = 1,
                                          add_noise = TRUE,
                                          seed = NULL) {
  stopifnot(nrow(biomarkers) > 0L)
  needed_cols <- c("rt_loc", "cv_loc", "sigma_rt", "sigma_cv", "intensity",
                   "case_prevalence", "control_prevalence",
                   "case_intensity_mult", "control_intensity_mult")
  missing <- setdiff(needed_cols, names(biomarkers))
  if (length(missing)) {
    stop("biomarkers is missing columns: ",
         paste(missing, collapse = ", "))
  }
  if (!is.null(seed)) set.seed(seed)

  N <- n_cases + n_controls
  y <- c(rep(2L, n_cases), rep(1L, n_controls))
  sample_names <- c(sprintf("case_%03d", seq_len(n_cases)),
                    sprintf("control_%03d", seq_len(n_controls)))
  names(y) <- sample_names

  samples <- vector("list", N)
  names(samples) <- sample_names

  n_biomarkers <- nrow(biomarkers)
  firing <- matrix(FALSE, nrow = N, ncol = n_biomarkers,
                   dimnames = list(sample_names, NULL))

  rows <- seq_len(H); cols <- seq_len(W)

  for (i in seq_len(N)) {
    is_case <- y[i] == 2L

    # --- Background: draw clean from background peak_params, no noise ---
    bg_pair <- generate_one_synthetic(
      background_peak_params, H = H, W = W,
      add_noise          = FALSE,
      size_jitter        = size_jitter,
      location_mode      = "empirical",
      location_jitter_rt = location_jitter_rt,
      location_jitter_cv = location_jitter_cv,
      attribute_mode     = "joint"
    )
    clean <- bg_pair$clean

    # --- Biomarker injection: Bernoulli per biomarker, add if fires ---
    for (b in seq_len(n_biomarkers)) {
      prev <- if (is_case) biomarkers$case_prevalence[b]
              else         biomarkers$control_prevalence[b]
      if (stats::runif(1) >= prev) next
      firing[i, b] <- TRUE

      mult <- if (is_case) biomarkers$case_intensity_mult[b]
              else         biomarkers$control_intensity_mult[b]

      mu_rt <- biomarkers$rt_loc[b] * H +
                 stats::rnorm(1, 0, location_jitter_rt)
      mu_cv <- biomarkers$cv_loc[b] * W +
                 stats::rnorm(1, 0, location_jitter_cv)
      mu_rt <- max(1, min(H, mu_rt))
      mu_cv <- max(1, min(W, mu_cv))

      # Same size-jitter + amplitude modulation as generate_one_synthetic
      # so injected peaks are visually indistinguishable from background peaks
      # (except for their location and intensity).
      scale_factor <- exp(stats::rnorm(1, 0, size_jitter))
      sig_rt <- max(biomarkers$sigma_rt[b] * scale_factor, 0.8)
      sig_cv <- max(biomarkers$sigma_cv[b] * scale_factor, 0.8)
      base_amp <- biomarkers$intensity[b] * mult
      amp <- base_amp * (0.5 + 0.5 * scale_factor)

      g_rt <- exp(-0.5 * ((rows - mu_rt) / sig_rt)^2)
      g_cv <- exp(-0.5 * ((cols - mu_cv) / sig_cv)^2)
      clean <- clean + amp * outer(g_rt, g_cv)
    }

    # --- Noise ---
    if (add_noise &&
        is.finite(background_peak_params$noise$mean) &&
        is.finite(background_peak_params$noise$sd)) {
      noise_sd <- max(1e-8, background_peak_params$noise$sd)
      noise <- matrix(
        abs(stats::rnorm(H * W,
                         mean = background_peak_params$noise$mean,
                         sd   = noise_sd)),
        nrow = H, ncol = W
      )
      Z <- clean + noise
    } else {
      Z <- clean
    }

    samples[[i]] <- list(
      path = paste0("simulated/", sample_names[i]),
      time = seq_len(H),   # placeholder RT axis
      cv   = seq_len(W),   # placeholder CV axis
      Z    = Z
    )
  }

  list(
    samples           = samples,
    y                 = y,
    ground_truth      = biomarkers,
    per_sample_firing = firing
  )
}
