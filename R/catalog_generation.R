# classydms: catalog-based synthetic generation and case/control simulation
# Uses a peak_catalog (from build_peak_catalog) rather than per-sample
# independent draws. This preserves real-data shared-compound structure:
# every synthetic sample draws from the same catalog of compound
# positions, with per-compound prevalence controlling how often each
# compound fires.

#' Generate one synthetic sample from a peak catalog
#'
#' For each compound in the catalog, a `Bernoulli(prevalence)` draw
#' decides whether the compound "fires" in this sample. If it fires, a
#' 2-D Gaussian peak is placed at `compound$rt_loc, compound$cv_loc`
#' with morphology bootstrapped from the compound's observed
#' `(sigma_rt, sigma_cv, intensity)` triples. Per-peak location and
#' size jitter apply on top. Folded-Normal noise is added at the end
#' when `add_noise = TRUE`.
#'
#' Compared to [generate_one_synthetic()], the catalog version reflects
#' the real-data property that the same compound positions appear
#' across many samples. Which specific compounds appear per sample is
#' random (via Bernoulli), but the compound *identities* and positions
#' are shared across all synthetic samples.
#'
#' @param catalog A `peak_catalog` object from [build_peak_catalog()].
#' @param H,W Output image dimensions in pixels.
#' @param add_noise Whether to add folded-Normal background noise.
#' @param size_jitter Per-peak log-Normal size scale factor SD.
#' @param location_jitter_rt,location_jitter_cv Per-peak location
#'   jitter SD (pixels).
#' @param noise_scale Multiplier on `catalog$noise$sd`. Default `1.0`.
#' @param prevalence_override Optional named numeric vector mapping
#'   `compound_id` (as character) to a replacement prevalence. Not
#'   present compounds use the catalog's own prevalence.
#' @param intensity_mult Optional named numeric vector mapping
#'   `compound_id` (as character) to an amplitude multiplier applied
#'   when the compound fires. Defaults to `1.0` for unlisted compounds.
#' @param compound_ids Optional integer vector of compound_ids to
#'   include; `NULL` (default) uses every catalog compound.
#'
#' @return List with `clean` (H x W), `noisy` (H x W, `= clean` when
#'   `add_noise = FALSE`), and `fired` (logical vector of length
#'   `nrow(catalog$compounds)` indicating which compounds actually
#'   appeared in this sample).
#' @export
generate_one_synthetic_from_catalog <- function(catalog, H, W,
                                                 add_noise = TRUE,
                                                 size_jitter = 0.15,
                                                 location_jitter_rt = 2,
                                                 location_jitter_cv = 1,
                                                 noise_scale = 1.0,
                                                 prevalence_override = NULL,
                                                 intensity_mult = NULL,
                                                 compound_ids = NULL) {
  stopifnot(inherits(catalog, "peak_catalog"))
  compounds <- catalog$compounds
  if (!is.null(compound_ids)) {
    compounds <- compounds[compounds$compound_id %in% compound_ids, ]
  }
  if (nrow(compounds) == 0L) {
    stop("generate_one_synthetic_from_catalog: no compounds to generate.")
  }
  n_compounds <- nrow(compounds)

  clean  <- matrix(0, nrow = H, ncol = W)
  rows   <- seq_len(H); cols <- seq_len(W)
  fired  <- logical(n_compounds)

  # Precompute lookup tables for override vectors
  po_names <- if (!is.null(prevalence_override))
                names(prevalence_override) else character(0)
  im_names <- if (!is.null(intensity_mult))
                names(intensity_mult) else character(0)

  for (i in seq_len(n_compounds)) {
    cid <- as.character(compounds$compound_id[i])
    prev <- if (cid %in% po_names) prevalence_override[[cid]]
            else                   compounds$prevalence[i]
    if (stats::runif(1) >= prev) next
    fired[i] <- TRUE

    mult <- if (cid %in% im_names) intensity_mult[[cid]] else 1.0

    # Bootstrap morphology triple from this compound's observations
    obs_intensity <- compounds$intensity_obs[[i]]
    n_obs <- length(obs_intensity)
    j <- sample.int(n_obs, size = 1L)
    sig_rt   <- compounds$sigma_rt_obs[[i]][j]
    sig_cv   <- compounds$sigma_cv_obs[[i]][j]
    base_amp <- obs_intensity[j]

    # Place with location jitter
    mu_rt <- compounds$rt_loc[i] * H +
               stats::rnorm(1, 0, location_jitter_rt)
    mu_cv <- compounds$cv_loc[i] * W +
               stats::rnorm(1, 0, location_jitter_cv)
    mu_rt <- max(1, min(H, mu_rt))
    mu_cv <- max(1, min(W, mu_cv))

    # Size jitter and amplitude modulation (same convention as
    # generate_one_synthetic so morphology looks consistent)
    scale_factor <- exp(stats::rnorm(1, 0, size_jitter))
    sig_rt <- max(sig_rt * scale_factor, 0.8)
    sig_cv <- max(sig_cv * scale_factor, 0.8)
    amp <- base_amp * mult * (0.5 + 0.5 * scale_factor)

    g_rt <- exp(-0.5 * ((rows - mu_rt) / sig_rt)^2)
    g_cv <- exp(-0.5 * ((cols - mu_cv) / sig_cv)^2)
    clean <- clean + amp * outer(g_rt, g_cv)
  }

  if (!add_noise || is.null(catalog$noise) ||
      !is.finite(catalog$noise$mean) || !is.finite(catalog$noise$sd)) {
    return(list(clean = clean, noisy = clean, fired = fired))
  }
  noise_sd <- max(1e-8, catalog$noise$sd * noise_scale)
  noise <- matrix(abs(stats::rnorm(H * W, mean = catalog$noise$mean,
                                    sd = noise_sd)),
                   nrow = H, ncol = W)
  list(clean = clean, noisy = clean + noise, fired = fired)
}

#' Pick a set of catalog compounds as biomarkers
#'
#' Selects `n_biomarkers` compound_ids uniformly at random from the
#' catalog, optionally filtered by prevalence range or intensity
#' quantile. Analogous to [holdout_biomarkers()] for the peak_params
#' flow.
#'
#' Unlike the peak_params version, no "protect radius" is needed —
#' catalog compounds already sit at discrete positions with no
#' overlap between them, so background compounds (all the non-selected
#' entries) can't land in biomarker regions.
#'
#' @param catalog A `peak_catalog` object.
#' @param n_biomarkers Number of biomarker compound_ids to pick.
#' @param min_prevalence,max_prevalence Only consider compounds whose
#'   catalog prevalence falls in this range. Defaults `[0, 1]`
#'   (any). For "rare-in-cohort" biomarkers set e.g. `max = 0.3`; for
#'   "common-in-cohort" set `min = 0.7`.
#' @param min_intensity_quantile Optional threshold on the median of
#'   each compound's `intensity_obs`. Same idea as in
#'   [holdout_biomarkers()].
#' @param seed Optional RNG seed.
#'
#' @return Integer vector of length `n_biomarkers` — the selected
#'   compound_ids.
#' @export
pick_catalog_biomarkers <- function(catalog, n_biomarkers,
                                     min_prevalence = 0,
                                     max_prevalence = 1,
                                     min_intensity_quantile = NULL,
                                     seed = NULL) {
  stopifnot(inherits(catalog, "peak_catalog"))
  if (!is.null(seed)) set.seed(seed)

  df <- catalog$compounds
  keep <- df$prevalence >= min_prevalence &
          df$prevalence <= max_prevalence
  if (!is.null(min_intensity_quantile)) {
    median_intensity <- vapply(df$intensity_obs, stats::median,
                                numeric(1))
    thr <- as.numeric(stats::quantile(median_intensity,
                                        min_intensity_quantile))
    keep <- keep & (median_intensity >= thr)
  }
  eligible <- df$compound_id[keep]
  if (length(eligible) < n_biomarkers) {
    stop("pick_catalog_biomarkers: only ", length(eligible),
         " compounds match the filters, but ", n_biomarkers,
         " biomarkers requested.")
  }
  sort(sample(eligible, size = n_biomarkers, replace = FALSE))
}

#' Build a case/control biomarker spec for catalog-based simulation
#'
#' Wraps selected compound_ids with the per-group prevalence and
#' intensity overrides used by [simulate_case_control_from_catalog()].
#'
#' @param compound_ids Integer vector of compound_ids designated as
#'   biomarkers (see [pick_catalog_biomarkers()]).
#' @param case_prevalence,control_prevalence Bernoulli probability of
#'   the biomarker firing in a case / control sample. Scalar
#'   (broadcast) or length-matched vector.
#' @param case_intensity_mult,control_intensity_mult Amplitude
#'   multiplier applied to the compound's bootstrapped intensity in
#'   case / control samples. Scalar or length-matched vector.
#'
#' @return A tibble with columns `compound_id`, `case_prevalence`,
#'   `control_prevalence`, `case_intensity_mult`,
#'   `control_intensity_mult`. Pass to
#'   [simulate_case_control_from_catalog()].
#' @export
build_biomarker_spec_catalog <- function(compound_ids,
                                          case_prevalence,
                                          control_prevalence,
                                          case_intensity_mult    = 1,
                                          control_intensity_mult = 1) {
  n <- length(compound_ids)
  bcast <- function(x, nm) {
    if (length(x) == 1L) return(rep(x, n))
    if (length(x) == n)  return(x)
    stop(nm, " must have length 1 or ", n, ", got ", length(x), ".")
  }
  tibble::tibble(
    compound_id            = as.integer(compound_ids),
    case_prevalence        = bcast(case_prevalence, "case_prevalence"),
    control_prevalence     = bcast(control_prevalence,
                                     "control_prevalence"),
    case_intensity_mult    = bcast(case_intensity_mult,
                                     "case_intensity_mult"),
    control_intensity_mult = bcast(control_intensity_mult,
                                     "control_intensity_mult")
  )
}

#' Simulate a case / control cohort from a peak catalog
#'
#' For each sample: determine group (case or control), then generate
#' via [generate_one_synthetic_from_catalog()] with the biomarker
#' spec's per-group prevalence and intensity overrides applied.
#' Non-biomarker compounds use the catalog's own prevalence in both
#' groups (background is shared). This directly matches the real-data
#' property that the *same* compound positions appear across most
#' samples, only differing in whether individual compounds fired.
#'
#' Compared to [simulate_case_control_cohort()] (the peak_params flow):
#' - Shared background: every sample draws from the same catalog
#'   compound positions rather than random-per-sample peak positions.
#' - No protect_radius: catalog compounds are already at discrete
#'   positions with no overlap between them.
#' - Biomarker specification is by `compound_id`, not by random draw
#'   from a pool.
#'
#' @param catalog A `peak_catalog` object.
#' @param biomarkers A tibble from [build_biomarker_spec_catalog()].
#' @param n_cases,n_controls Number of samples per group.
#' @param H,W Output image dimensions in pixels.
#' @param size_jitter,location_jitter_rt,location_jitter_cv Per-peak
#'   jitter parameters, forwarded to the underlying generator.
#' @param add_noise Whether to add background noise.
#' @param noise_scale Multiplier on catalog noise SD.
#' @param seed Optional RNG seed.
#'
#' @return List with `samples`, `y`, `ground_truth` (the biomarker
#'   spec), and `per_sample_firing` (logical matrix of `n_samples` x
#'   `n_biomarkers`; whether each biomarker fired in each sample).
#'   `samples` matches [load_dms_directory()]'s shape so downstream
#'   pipeline steps work unchanged.
#' @export
simulate_case_control_from_catalog <- function(catalog,
                                                 biomarkers,
                                                 n_cases, n_controls,
                                                 H, W,
                                                 size_jitter = 0.15,
                                                 location_jitter_rt = 2,
                                                 location_jitter_cv = 1,
                                                 add_noise = TRUE,
                                                 noise_scale = 1.0,
                                                 seed = NULL) {
  stopifnot(inherits(catalog, "peak_catalog"))
  needed <- c("compound_id", "case_prevalence", "control_prevalence",
              "case_intensity_mult", "control_intensity_mult")
  missing <- setdiff(needed, names(biomarkers))
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

  # Precompute the per-group overrides once
  cid_chr <- as.character(biomarkers$compound_id)
  case_prev_override <- stats::setNames(
    biomarkers$case_prevalence, cid_chr)
  ctrl_prev_override <- stats::setNames(
    biomarkers$control_prevalence, cid_chr)
  case_int_mult <- stats::setNames(
    biomarkers$case_intensity_mult, cid_chr)
  ctrl_int_mult <- stats::setNames(
    biomarkers$control_intensity_mult, cid_chr)

  samples <- vector("list", N)
  names(samples) <- sample_names
  n_biomarkers <- nrow(biomarkers)
  firing <- matrix(FALSE, nrow = N, ncol = n_biomarkers,
                   dimnames = list(sample_names, NULL))

  # Which rows of catalog$compounds correspond to biomarker ids?
  biom_rows <- match(biomarkers$compound_id,
                      catalog$compounds$compound_id)
  if (any(is.na(biom_rows))) {
    bad <- biomarkers$compound_id[is.na(biom_rows)]
    stop("biomarkers references compound_id(s) not in the catalog: ",
         paste(bad, collapse = ", "))
  }

  for (i in seq_len(N)) {
    is_case <- y[i] == 2L
    prev_override <- if (is_case) case_prev_override
                     else         ctrl_prev_override
    int_mult      <- if (is_case) case_int_mult
                     else         ctrl_int_mult

    result <- generate_one_synthetic_from_catalog(
      catalog, H = H, W = W,
      add_noise           = add_noise,
      size_jitter         = size_jitter,
      location_jitter_rt  = location_jitter_rt,
      location_jitter_cv  = location_jitter_cv,
      noise_scale         = noise_scale,
      prevalence_override = prev_override,
      intensity_mult      = int_mult
    )

    # Record biomarker firing status for this sample
    firing[i, ] <- result$fired[biom_rows]

    samples[[i]] <- list(
      path = paste0("simulated/", sample_names[i]),
      time = seq_len(H),
      cv   = seq_len(W),
      Z    = result$noisy
    )
  }

  list(
    samples           = samples,
    y                 = y,
    ground_truth      = biomarkers,
    per_sample_firing = firing
  )
}

#' Sanity-max scenario for catalog-based case/control simulation
#'
#' Analogous to [biomarkers_sanity_max()] for the peak_params flow.
#' Picks 5 catalog compounds whose intensity is in the top 5% of
#' compound median intensities and whose catalog prevalence is at
#' least `min_base_prevalence` (default 0.5). Sets case prevalence to
#' 1.0, control prevalence to 0.0, and case_intensity_mult to 10x.
#' Designed to test that the full catalog-based pipeline can recover
#' an obviously injected signal.
#'
#' @param catalog A `peak_catalog` object.
#' @param min_base_prevalence Only pick biomarkers whose catalog
#'   prevalence is >= this. Default `0.5`.
#' @param seed Optional RNG seed.
#' @return A biomarker spec tibble ready for
#'   [simulate_case_control_from_catalog()].
#' @export
catalog_biomarkers_sanity_max <- function(catalog,
                                            min_base_prevalence = 0.5,
                                            seed = NULL) {
  compound_ids <- pick_catalog_biomarkers(
    catalog,
    n_biomarkers = 5L,
    min_prevalence = min_base_prevalence,
    max_prevalence = 1.0,
    min_intensity_quantile = 0.95,
    seed = seed
  )
  build_biomarker_spec_catalog(
    compound_ids,
    case_prevalence        = 1.0,
    control_prevalence     = 0.0,
    case_intensity_mult    = 10.0,
    control_intensity_mult = 1.0
  )
}
