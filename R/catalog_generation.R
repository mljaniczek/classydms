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

  # Placement uses FRACTION coordinates so peaks land at fraction * H / W.
  # For physical-mode catalogs, rt_loc / cv_loc are in seconds / volts and
  # rt_frac / cv_frac carry the fraction. Read coord_mode to pick the
  # right columns. Older catalogs (no coord_mode field) are treated as
  # fraction for backward compatibility.
  coord_mode <- catalog$parameters$coord_mode %||% "fraction"
  if (coord_mode == "physical") {
    if (is.null(compounds$rt_frac) || is.null(compounds$cv_frac)) {
      stop("generate_one_synthetic_from_catalog: physical-mode catalog ",
           "must carry rt_frac and cv_frac columns. Rebuild the catalog ",
           "with a recent build_peak_catalog().")
    }
    rt_place <- compounds$rt_frac
    cv_place <- compounds$cv_frac
  } else {
    rt_place <- compounds$rt_loc
    cv_place <- compounds$cv_loc
  }

  clean  <- matrix(0, nrow = H, ncol = W)
  rows   <- seq_len(H); cols <- seq_len(W)
  fired  <- logical(n_compounds)
  # Track per-fired-compound placement so downstream code can subtract or
  # add peaks by exact morphology (used for cofire post-processing).
  placements_cid   <- integer(0)
  placements_mu_rt <- numeric(0)
  placements_mu_cv <- numeric(0)
  placements_sig_rt <- numeric(0)
  placements_sig_cv <- numeric(0)
  placements_amp   <- numeric(0)

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

    # If the compound's placement fraction is non-finite (a bad catalog
    # entry that slipped through), skip it. This can happen if the
    # catalog was built from peak_params containing NAs in
    # rt_loc_raw / cv_loc_raw and the median aggregation produced NaN.
    if (!is.finite(rt_place[i]) || !is.finite(cv_place[i])) next

    mult <- if (cid %in% im_names) intensity_mult[[cid]] else 1.0

    # Bootstrap morphology triple from this compound's observations.
    # Skip the compound (don't count as fired) if we can't find a valid
    # finite triple in a bounded number of tries — the alternative is a
    # NaN pixel that will crash normalization downstream.
    obs_intensity <- compounds$intensity_obs[[i]]
    n_obs <- length(obs_intensity)
    if (n_obs == 0L) next
    got_valid <- FALSE
    for (draw_attempt in seq_len(5L)) {
      j <- sample.int(n_obs, size = 1L)
      sig_rt   <- compounds$sigma_rt_obs[[i]][j]
      sig_cv   <- compounds$sigma_cv_obs[[i]][j]
      base_amp <- obs_intensity[j]
      if (is.finite(sig_rt) && sig_rt > 0 &&
          is.finite(sig_cv) && sig_cv > 0 &&
          is.finite(base_amp) && base_amp > 0) {
        got_valid <- TRUE; break
      }
    }
    if (!got_valid) next
    fired[i] <- TRUE

    # Place with location jitter (rt_place / cv_place are always in
    # fraction — see coord_mode branch above).
    mu_rt <- rt_place[i] * H +
               stats::rnorm(1, 0, location_jitter_rt)
    mu_cv <- cv_place[i] * W +
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
    placements_cid    <- c(placements_cid,    compounds$compound_id[i])
    placements_mu_rt  <- c(placements_mu_rt,  mu_rt)
    placements_mu_cv  <- c(placements_mu_cv,  mu_cv)
    placements_sig_rt <- c(placements_sig_rt, sig_rt)
    placements_sig_cv <- c(placements_sig_cv, sig_cv)
    placements_amp    <- c(placements_amp,    amp)
  }

  placements <- tibble::tibble(
    compound_id = placements_cid,
    mu_rt = placements_mu_rt, mu_cv = placements_mu_cv,
    sig_rt = placements_sig_rt, sig_cv = placements_sig_cv,
    amp = placements_amp
  )

  if (!add_noise || is.null(catalog$noise) ||
      !is.finite(catalog$noise$mean) || !is.finite(catalog$noise$sd)) {
    stopifnot(
      "generate_one_synthetic_from_catalog produced non-finite clean pixels" =
        all(is.finite(clean)))
    return(list(clean = clean, noisy = clean, fired = fired,
                placements = placements))
  }
  noise_sd   <- max(1e-8, catalog$noise$sd * noise_scale)
  noise_mean <- catalog$noise$mean
  # rnorm silently returns NaN if mean or sd is non-finite; guard here.
  if (!is.finite(noise_mean)) noise_mean <- 0
  if (!is.finite(noise_sd))   noise_sd   <- 1e-8
  noise <- matrix(abs(stats::rnorm(H * W, mean = noise_mean, sd = noise_sd)),
                   nrow = H, ncol = W)
  noisy <- clean + noise
  stopifnot(
    "generate_one_synthetic_from_catalog produced non-finite pixels" =
      all(is.finite(clean)) && all(is.finite(noisy)))
  list(clean = clean, noisy = noisy, fired = fired,
       placements = placements)
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
    # Length-safe + NA-safe median. A compound with an empty intensity_obs
    # (older catalogs pre-fix, or catalogs loaded from disk) would give
    # median() a numeric(0) and yield NaN — filter these out of the
    # eligibility set rather than letting them poison quantile().
    median_intensity <- vapply(df$intensity_obs, function(x) {
      x <- x[is.finite(x) & x > 0]
      if (length(x) == 0L) NA_real_ else stats::median(x)
    }, numeric(1))
    finite_mask <- is.finite(median_intensity)
    thr <- as.numeric(stats::quantile(median_intensity[finite_mask],
                                        min_intensity_quantile,
                                        na.rm = TRUE))
    keep <- keep & finite_mask & (median_intensity >= thr)
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
#' @param cofire Optional list of dependency groups modelling co-firing
#'   structure between compounds. Each element is a list with EITHER
#'   `anchor` (a single `compound_id`), `dependents` (integer vector of
#'   `compound_id`s), and `p_dependent_given_anchor` (scalar in `[0,1]`
#'   — probability each dependent fires when the anchor does), OR
#'   `members` (integer vector of `compound_id`s) and `p_joint` (scalar
#'   in `[0,1]` — probability the whole clique fires together). Applied
#'   as post-processing on top of the independent Bernoulli firing:
#'   - Dependency group: if anchor did not fire, every dependent that
#'     independently fired is REMOVED (its Gaussian subtracted). If the
#'     anchor fired, every dependent that did not fire independently gets
#'     an additional Bernoulli(`p_dependent_given_anchor`) roll; on
#'     success the compound fires and its Gaussian is added.
#'   - Symmetric clique: one joint Bernoulli(`p_joint`) roll decides the
#'     whole group. On success every non-fired member is added; on
#'     failure every fired member is removed.
#'   The final `per_sample_firing` matrix reflects the post-cofire state.
#' @param off_catalog_peaks Optional list controlling injection of random
#'   peaks NOT present in the catalog. Fields:
#'   - `n_cases`, `n_controls`: integer peaks to inject per case /
#'     control sample.
#'   - `rt_range`, `cv_range`: length-2 numerics in [0,1] giving the
#'     fractional (RT, CV) placement region.
#'   - `intensity_range`: length-2 numeric giving uniform bounds for
#'     the injected peak amplitude.
#'   - `exclusion_eps_rt`, `exclusion_eps_cv`: fractional radius; a
#'     candidate position is rejected if within this box of any catalog
#'     compound. Defaults `(0.02, 0.05)`.
#'   Sigmas are drawn from the pooled catalog sigma distribution so peak
#'   shapes match the catalog's morphology. Used for stress-testing an
#'   encoder / classifier on peaks the pretraining never saw. Peak
#'   positions are recorded on the returned `off_catalog_positions`
#'   list (one entry per sample).
#' @param seed Optional RNG seed.
#'
#' @return List with `samples`, `y`, `ground_truth` (the biomarker
#'   spec), `per_sample_firing` (logical matrix of `n_samples` x
#'   `n_biomarkers`; whether each biomarker fired in each sample AFTER
#'   any cofire post-processing), `cofire` (the cofire spec that was
#'   applied, or `NULL`), and `off_catalog_positions` (list of length
#'   `n_samples` recording injected positions, or `NULL`).
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
                                                 cofire = NULL,
                                                 off_catalog_peaks = NULL,
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

  # ---- cofire validation ----
  if (!is.null(cofire)) {
    if (!is.list(cofire) || !all(vapply(cofire, is.list, logical(1)))) {
      stop("cofire must be a list of lists (one per group).")
    }
    for (g_idx in seq_along(cofire)) {
      g <- cofire[[g_idx]]
      is_dep_group <- !is.null(g$anchor)
      is_clique    <- !is.null(g$members)
      if (is_dep_group + is_clique != 1L) {
        stop("cofire[[", g_idx, "]] must have EITHER anchor+dependents ",
             "OR members (not both, not neither).")
      }
      if (is_dep_group) {
        if (is.null(g$dependents) || is.null(g$p_dependent_given_anchor))
          stop("cofire[[", g_idx,
               "]] dependency group requires anchor, dependents, ",
               "p_dependent_given_anchor.")
        cids <- c(g$anchor, g$dependents)
      } else {
        if (is.null(g$p_joint))
          stop("cofire[[", g_idx, "]] clique requires members, p_joint.")
        cids <- g$members
      }
      if (!all(cids %in% catalog$compounds$compound_id)) {
        bad <- setdiff(cids, catalog$compounds$compound_id)
        stop("cofire[[", g_idx,
             "]] references compound_id(s) not in the catalog: ",
             paste(bad, collapse = ", "))
      }
    }
  }

  # ---- off_catalog_peaks validation + preparation ----
  if (!is.null(off_catalog_peaks)) {
    ocp <- off_catalog_peaks
    ocp$exclusion_eps_rt <- ocp$exclusion_eps_rt %||% 0.02
    ocp$exclusion_eps_cv <- ocp$exclusion_eps_cv %||% 0.05
    for (nm in c("n_cases", "n_controls", "rt_range", "cv_range",
                  "intensity_range")) {
      if (is.null(ocp[[nm]]))
        stop("off_catalog_peaks$", nm, " is required.")
    }
    # Pool sigma distributions. Filter to finite positive values so a
    # catalog with stale NAs in the list-columns (pre-fix build or
    # externally-provided catalog) can't seed off-catalog peaks with
    # NA sigma that render as NaN Gaussians.
    srp <- unlist(catalog$compounds$sigma_rt_obs)
    scp <- unlist(catalog$compounds$sigma_cv_obs)
    ocp$sigma_rt_pool <- srp[is.finite(srp) & srp > 0]
    ocp$sigma_cv_pool <- scp[is.finite(scp) & scp > 0]
    if (length(ocp$sigma_rt_pool) == 0L || length(ocp$sigma_cv_pool) == 0L)
      stop("off_catalog_peaks: catalog has no usable sigma observations ",
           "to draw peak shapes from.")
    # Catalog compound positions in fraction space (for exclusion)
    cmp <- catalog$compounds
    ocp$cat_rt_frac <- if (!is.null(cmp$rt_frac)) cmp$rt_frac else cmp$rt_loc
    ocp$cat_cv_frac <- if (!is.null(cmp$cv_frac)) cmp$cv_frac else cmp$cv_loc
    # Filter to finite so exclusion check doesn't compare against NA
    keep_pos <- is.finite(ocp$cat_rt_frac) & is.finite(ocp$cat_cv_frac)
    ocp$cat_rt_frac <- ocp$cat_rt_frac[keep_pos]
    ocp$cat_cv_frac <- ocp$cat_cv_frac[keep_pos]
  } else {
    ocp <- NULL
  }

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

  # Coord-mode-aware placement columns (same convention as
  # generate_one_synthetic_from_catalog).
  coord_mode <- catalog$parameters$coord_mode %||% "fraction"
  cmp <- catalog$compounds
  rt_place <- if (coord_mode == "physical") cmp$rt_frac else cmp$rt_loc
  cv_place <- if (coord_mode == "physical") cmp$cv_frac else cmp$cv_loc

  off_catalog_positions <- if (!is.null(ocp)) vector("list", N) else NULL
  if (!is.null(off_catalog_positions))
    names(off_catalog_positions) <- sample_names

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

    Z <- result$clean          # work on clean, add noise last
    fired <- result$fired
    placements <- result$placements

    # ---- cofire post-processing ----
    if (!is.null(cofire)) {
      cid_to_row <- stats::setNames(seq_len(nrow(catalog$compounds)),
                                     as.character(catalog$compounds$compound_id))
      # Fast index into placements by compound_id (may hold >1 hit if
      # generator was called with duplicates; here always 1 or 0).
      placements_by_cid <- split(seq_len(nrow(placements)),
                                   placements$compound_id)

      subtract_compound <- function(cid) {
        idxs <- placements_by_cid[[as.character(cid)]]
        if (is.null(idxs)) return()
        for (p in idxs) {
          g_rt <- exp(-0.5 * ((seq_len(H) - placements$mu_rt[p]) /
                                placements$sig_rt[p])^2)
          g_cv <- exp(-0.5 * ((seq_len(W) - placements$mu_cv[p]) /
                                placements$sig_cv[p])^2)
          Z <<- Z - placements$amp[p] * outer(g_rt, g_cv)
        }
      }
      add_compound <- function(cid) {
        crow <- cid_to_row[[as.character(cid)]]
        if (is.na(crow)) return()
        # Skip if placement fraction non-finite (catalog defense).
        if (!is.finite(rt_place[crow]) || !is.finite(cv_place[crow])) return()
        # Draw fresh morphology (same convention as generator, with the
        # same resample-on-bad-triple loop so a stale NA in the per-obs
        # lists doesn't render a NaN Gaussian and break subtract later).
        obs_intensity <- catalog$compounds$intensity_obs[[crow]]
        n_obs <- length(obs_intensity)
        if (n_obs == 0L) return()
        got_valid <- FALSE
        for (draw_attempt in seq_len(5L)) {
          j <- sample.int(n_obs, size = 1L)
          sig_rt <- catalog$compounds$sigma_rt_obs[[crow]][j]
          sig_cv <- catalog$compounds$sigma_cv_obs[[crow]][j]
          base_amp <- obs_intensity[j]
          if (is.finite(sig_rt) && sig_rt > 0 &&
              is.finite(sig_cv) && sig_cv > 0 &&
              is.finite(base_amp) && base_amp > 0) {
            got_valid <- TRUE; break
          }
        }
        if (!got_valid) return()
        mult <- if (as.character(cid) %in% names(int_mult))
                  int_mult[[as.character(cid)]] else 1.0
        mu_rt <- rt_place[crow] * H + stats::rnorm(1, 0, location_jitter_rt)
        mu_cv <- cv_place[crow] * W + stats::rnorm(1, 0, location_jitter_cv)
        mu_rt <- max(1, min(H, mu_rt))
        mu_cv <- max(1, min(W, mu_cv))
        scale_factor <- exp(stats::rnorm(1, 0, size_jitter))
        sig_rt <- max(sig_rt * scale_factor, 0.8)
        sig_cv <- max(sig_cv * scale_factor, 0.8)
        amp <- base_amp * mult * (0.5 + 0.5 * scale_factor)
        g_rt <- exp(-0.5 * ((seq_len(H) - mu_rt) / sig_rt)^2)
        g_cv <- exp(-0.5 * ((seq_len(W) - mu_cv) / sig_cv)^2)
        Z <<- Z + amp * outer(g_rt, g_cv)
      }

      for (g in cofire) {
        if (!is.null(g$anchor)) {
          # Dependency group
          anchor_row <- cid_to_row[[as.character(g$anchor)]]
          anchor_fired <- fired[anchor_row]
          for (dcid in g$dependents) {
            drow <- cid_to_row[[as.character(dcid)]]
            if (is.na(drow)) next
            if (!anchor_fired && fired[drow]) {
              subtract_compound(dcid); fired[drow] <- FALSE
            } else if (anchor_fired && !fired[drow]) {
              if (stats::runif(1) < g$p_dependent_given_anchor) {
                add_compound(dcid); fired[drow] <- TRUE
              }
            }
          }
        } else {
          # Symmetric clique
          fire_all <- stats::runif(1) < g$p_joint
          for (mcid in g$members) {
            mrow <- cid_to_row[[as.character(mcid)]]
            if (is.na(mrow)) next
            if (fire_all && !fired[mrow]) {
              add_compound(mcid); fired[mrow] <- TRUE
            } else if (!fire_all && fired[mrow]) {
              subtract_compound(mcid); fired[mrow] <- FALSE
            }
          }
        }
      }
    }

    # ---- off-catalog injection ----
    if (!is.null(ocp)) {
      n_inj <- if (is_case) ocp$n_cases else ocp$n_controls
      if (n_inj > 0L) {
        pos_df <- inject_off_catalog(Z, ocp, H, W, n_inj,
                                       size_jitter,
                                       location_jitter_rt,
                                       location_jitter_cv)
        Z <- pos_df$Z
        off_catalog_positions[[i]] <- pos_df$positions
      } else {
        off_catalog_positions[[i]] <- data.frame(
          mu_rt = numeric(0), mu_cv = numeric(0),
          sig_rt = numeric(0), sig_cv = numeric(0), amp = numeric(0))
      }
    }

    # Re-apply noise last (subtract old noise, add fresh) so injections
    # don't leave stale noise unrelated to the modified peaks.
    if (add_noise && !is.null(catalog$noise) &&
        is.finite(catalog$noise$mean) && is.finite(catalog$noise$sd)) {
      noise_sd   <- max(1e-8, catalog$noise$sd * noise_scale)
      noise_mean <- catalog$noise$mean
      if (!is.finite(noise_sd))   noise_sd   <- 1e-8
      if (!is.finite(noise_mean)) noise_mean <- 0
      noise <- matrix(abs(stats::rnorm(H * W,
                                        mean = noise_mean,
                                        sd = noise_sd)),
                       nrow = H, ncol = W)
      Z <- Z + noise
    }
    stopifnot(
      "simulate_case_control_from_catalog produced non-finite pixels" =
        all(is.finite(Z)))

    firing[i, ] <- fired[biom_rows]

    samples[[i]] <- list(
      path = paste0("simulated/", sample_names[i]),
      time = seq_len(H),
      cv   = seq_len(W),
      Z    = Z
    )
  }

  list(
    samples               = samples,
    y                     = y,
    ground_truth          = biomarkers,
    per_sample_firing     = firing,
    cofire                = cofire,
    off_catalog_positions = off_catalog_positions
  )
}

# Internal helper: inject n random off-catalog peaks and return updated Z
# plus a positions data frame.
inject_off_catalog <- function(Z, ocp, H, W, n,
                                 size_jitter,
                                 location_jitter_rt,
                                 location_jitter_cv) {
  max_attempts <- 50L
  positions <- data.frame(mu_rt = numeric(0), mu_cv = numeric(0),
                           sig_rt = numeric(0), sig_cv = numeric(0),
                           amp = numeric(0))
  for (k in seq_len(n)) {
    accepted <- FALSE
    for (attempt in seq_len(max_attempts)) {
      rt_frac <- stats::runif(1, ocp$rt_range[1], ocp$rt_range[2])
      cv_frac <- stats::runif(1, ocp$cv_range[1], ocp$cv_range[2])
      # Reject if within exclusion box of any catalog compound
      within_rt <- abs(ocp$cat_rt_frac - rt_frac) < ocp$exclusion_eps_rt
      within_cv <- abs(ocp$cat_cv_frac - cv_frac) < ocp$exclusion_eps_cv
      if (!any(within_rt & within_cv)) {
        accepted <- TRUE; break
      }
    }
    if (!accepted) next   # give up on this peak; sparse catalog corners rare
    sig_rt <- sample(ocp$sigma_rt_pool, size = 1L)
    sig_cv <- sample(ocp$sigma_cv_pool, size = 1L)
    scale_factor <- exp(stats::rnorm(1, 0, size_jitter))
    sig_rt <- max(sig_rt * scale_factor, 0.8)
    sig_cv <- max(sig_cv * scale_factor, 0.8)
    mu_rt <- rt_frac * H + stats::rnorm(1, 0, location_jitter_rt)
    mu_cv <- cv_frac * W + stats::rnorm(1, 0, location_jitter_cv)
    mu_rt <- max(1, min(H, mu_rt))
    mu_cv <- max(1, min(W, mu_cv))
    amp <- stats::runif(1, ocp$intensity_range[1], ocp$intensity_range[2])
    g_rt <- exp(-0.5 * ((seq_len(H) - mu_rt) / sig_rt)^2)
    g_cv <- exp(-0.5 * ((seq_len(W) - mu_cv) / sig_cv)^2)
    Z <- Z + amp * outer(g_rt, g_cv)
    positions <- rbind(positions,
      data.frame(mu_rt = mu_rt, mu_cv = mu_cv,
                  sig_rt = sig_rt, sig_cv = sig_cv, amp = amp))
  }
  list(Z = Z, positions = positions)
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

#' Universal-backdrop biomarker preset for catalog simulation
#'
#' Two-tier sample composition. Anchor layer: catalog compounds with
#' prevalence at or above `anchor_prevalence_threshold` (default 0.9)
#' are treated as universal-backdrop and appear in every simulated
#' sample at their catalog prevalence (typically ~1.0). Biomarker
#' layer: `n_biomarkers` non-anchor compounds are picked from the top
#' intensity quantile and given case/control Bernoulli prevalence
#' overrides, providing the classifier signal.
#'
#' Non-biomarker, non-anchor compounds keep their catalog empirical
#' prevalence (heterogeneous — some rare, some moderately common).
#'
#' Biomarkers are picked from BELOW the anchor threshold so the signal
#' stands out against the invariant backdrop. Filtering to top
#' `min_intensity_quantile` (default 0.95) ensures they are strong
#' enough to be reliably detected.
#'
#' @param catalog A `peak_catalog` object.
#' @param n_biomarkers Number of biomarker compounds. Default `5L`.
#' @param anchor_prevalence_threshold Compounds with catalog prevalence
#'   `>=` this are anchors and excluded from the biomarker pool.
#'   Default `0.9`. On a universal catalog (built via
#'   `build_universal_catalog()`) this typically selects the ~150
#'   compounds present across all sources.
#' @param biomarker_min_prevalence Minimum catalog prevalence for a
#'   non-anchor compound to be eligible as a biomarker. Default `0.05`
#'   — screens out ultra-rare compounds that may be jitter artifacts.
#' @param min_intensity_quantile Intensity quantile threshold applied
#'   to the compound's median observed intensity. Default `0.95`.
#' @param case_prevalence,control_prevalence Bernoulli firing
#'   probability of each biomarker in case / control samples.
#'   Defaults: `0.9` / `0.1`.
#' @param case_intensity_mult,control_intensity_mult Amplitude
#'   multiplier applied to biomarker peaks. Defaults: `1.2` / `1.0`.
#' @param seed Optional RNG seed.
#' @return A biomarker spec tibble (as
#'   [build_biomarker_spec_catalog()]) plus a `"n_anchors"` attribute
#'   recording how many anchor compounds were identified. Pass to
#'   [simulate_case_control_from_catalog()].
#' @export
catalog_biomarkers_universal_backdrop <- function(
  catalog,
  n_biomarkers = 5L,
  anchor_prevalence_threshold = 0.9,
  biomarker_min_prevalence = 0.05,
  min_intensity_quantile = 0.95,
  case_prevalence = 0.9,
  control_prevalence = 0.1,
  case_intensity_mult = 1.2,
  control_intensity_mult = 1.0,
  seed = NULL
) {
  stopifnot(inherits(catalog, "peak_catalog"))
  n_anchors <- sum(catalog$compounds$prevalence >=
                    anchor_prevalence_threshold)

  compound_ids <- pick_catalog_biomarkers(
    catalog,
    n_biomarkers = n_biomarkers,
    min_prevalence = biomarker_min_prevalence,
    # Strict inequality-lite: biomarkers must be BELOW the anchor
    # threshold so they aren't already always-on.
    max_prevalence = anchor_prevalence_threshold - 1e-9,
    min_intensity_quantile = min_intensity_quantile,
    seed = seed
  )
  spec <- build_biomarker_spec_catalog(
    compound_ids,
    case_prevalence        = case_prevalence,
    control_prevalence     = control_prevalence,
    case_intensity_mult    = case_intensity_mult,
    control_intensity_mult = control_intensity_mult
  )
  attr(spec, "n_anchors") <- n_anchors
  attr(spec, "anchor_prevalence_threshold") <- anchor_prevalence_threshold
  spec
}

#' Generate one three-layer synthetic sample for encoder pretraining
#'
#' Builds a synthetic GC-DMS image by stacking three layers of peaks:
#'
#' 1. **Anchor backbone.** Every compound in `anchor_ids` fires
#'    unconditionally, with its own catalog morphology and normal
#'    per-peak location / size jitter. These are the stable ~100
#'    highly-prevalent compounds you want the encoder to treat as the
#'    baseline landscape.
#' 2. **Variable catalog compounds** (optional). Per-sample random
#'    subset of the catalog's less-prevalent compounds — filtered by
#'    `variable_config$prevalence_range` and (optionally) excluding the
#'    anchor pool. Simulates diet-related, host-genetics-related, and
#'    other per-person volatiles. Intensity scaled by a per-peak factor
#'    drawn from `variable_config$intensity_scale`.
#' 3. **Contamination** (optional). Uniform-random peaks placed at
#'    positions NOT overlapping catalog compounds. Simulates sampling
#'    apparatus contamination, environmental exposure, and truly
#'    unfamiliar peaks. Sigmas drawn from the pooled catalog sigma
#'    distribution so shapes look realistic.
#'
#' @param catalog A `peak_catalog` object.
#' @param H,W Output image dimensions in pixels.
#' @param anchor_ids Integer vector of compound_ids to force-fire.
#'   Typically the ~100 anchors from the universal catalog with
#'   prevalence >= 0.9.
#' @param variable_config Optional list controlling the variable layer.
#'   Fields:
#'   - `n_per_sample`: scalar (fixed count) or length-2 numeric (uniform
#'     range) for how many variable compounds to draw. Default
#'     `c(30, 80)`.
#'   - `prevalence_range`: length-2 numeric giving the catalog
#'     prevalence range eligible for the variable pool. Default
#'     `c(0.05, 0.5)`.
#'   - `intensity_scale`: length-2 numeric giving amplitude multiplier
#'     bounds. Default `c(0.3, 1.0)` (variable compounds slightly
#'     dimmer than anchors on average).
#'   - `exclude_anchors`: `TRUE` (default) — remove anchor_ids from the
#'     variable pool so anchors don't get drawn twice.
#'   Pass `NULL` to skip this layer.
#' @param contamination_config Optional list controlling contamination
#'   injection. Same fields as `off_catalog_peaks` in
#'   [simulate_case_control_from_catalog()] except with `n_per_sample`
#'   instead of `n_cases` / `n_controls`. Default field values:
#'   `n_per_sample = c(5, 30)`, `rt_range = c(0.05, 0.95)`,
#'   `cv_range = c(0.05, 0.95)`, `intensity_range = c(0.1, 0.5)`,
#'   `exclusion_eps_rt = 0.02`, `exclusion_eps_cv = 0.05`. Pass `NULL`
#'   to skip this layer.
#' @param add_noise Whether to add folded-Normal background noise on
#'   top of the three layers.
#' @param size_jitter,location_jitter_rt,location_jitter_cv Per-peak
#'   jitter parameters (same convention as
#'   [generate_one_synthetic_from_catalog()]).
#' @param noise_scale Multiplier on catalog noise SD.
#' @return List with `clean` (H x W), `noisy` (H x W), and per-layer
#'   composition metadata (`anchor_ids`, `variable_ids`,
#'   `contamination_positions`).
#' @export
generate_noisy_pretrain_sample <- function(catalog, H, W, anchor_ids,
                                             variable_config = NULL,
                                             contamination_config = NULL,
                                             add_noise = TRUE,
                                             size_jitter = 0.15,
                                             location_jitter_rt = 2,
                                             location_jitter_cv = 1,
                                             noise_scale = 1.0) {
  stopifnot(inherits(catalog, "peak_catalog"))
  cmp <- catalog$compounds
  if (!all(anchor_ids %in% cmp$compound_id)) {
    bad <- setdiff(anchor_ids, cmp$compound_id)
    stop("anchor_ids references compound_id(s) not in the catalog: ",
         paste(head(bad, 5), collapse = ", "),
         if (length(bad) > 5) "..." else "")
  }

  # Layer 2 selection — draw variable compound IDs
  variable_ids <- integer(0)
  variable_int_mult <- setNames(numeric(0), character(0))
  if (!is.null(variable_config)) {
    vc <- variable_config
    vc$n_per_sample     <- vc$n_per_sample     %||% c(30, 80)
    vc$prevalence_range <- vc$prevalence_range %||% c(0.05, 0.5)
    vc$intensity_scale  <- vc$intensity_scale  %||% c(0.3, 1.0)
    vc$exclude_anchors  <- vc$exclude_anchors  %||% TRUE
    n_var <- if (length(vc$n_per_sample) == 1L) as.integer(vc$n_per_sample)
             else                                as.integer(stats::runif(1,
                                                    vc$n_per_sample[1],
                                                    vc$n_per_sample[2] + 1))
    eligible_mask <- cmp$prevalence >= vc$prevalence_range[1] &
                     cmp$prevalence <= vc$prevalence_range[2]
    if (vc$exclude_anchors)
      eligible_mask <- eligible_mask & !(cmp$compound_id %in% anchor_ids)
    eligible <- cmp$compound_id[eligible_mask]
    if (length(eligible) == 0L) {
      warning("generate_noisy_pretrain_sample: no compounds match ",
              "variable_config$prevalence_range; skipping variable layer.")
    } else {
      n_take <- min(n_var, length(eligible))
      variable_ids <- sample(eligible, size = n_take, replace = FALSE)
      variable_int_mult <- stats::setNames(
        stats::runif(length(variable_ids),
                      vc$intensity_scale[1], vc$intensity_scale[2]),
        as.character(variable_ids))
    }
  }

  fire_ids <- c(anchor_ids, variable_ids)
  # Force firing: prevalence_override = 1.0 for each. Anchor intensity
  # multiplier is 1.0; variable ones use their sampled scale.
  prev_ov <- stats::setNames(rep(1.0, length(fire_ids)),
                              as.character(fire_ids))
  int_ov  <- stats::setNames(rep(1.0, length(fire_ids)),
                              as.character(fire_ids))
  if (length(variable_int_mult))
    int_ov[names(variable_int_mult)] <- variable_int_mult

  base_res <- generate_one_synthetic_from_catalog(
    catalog, H = H, W = W,
    add_noise           = FALSE,
    size_jitter         = size_jitter,
    location_jitter_rt  = location_jitter_rt,
    location_jitter_cv  = location_jitter_cv,
    prevalence_override = prev_ov,
    intensity_mult      = int_ov,
    compound_ids        = fire_ids
  )
  Z <- base_res$clean

  # Layer 3 — contamination injection
  contamination_positions <- NULL
  if (!is.null(contamination_config)) {
    cc <- contamination_config
    cc$n_per_sample       <- cc$n_per_sample       %||% c(5, 30)
    cc$rt_range           <- cc$rt_range           %||% c(0.05, 0.95)
    cc$cv_range           <- cc$cv_range           %||% c(0.05, 0.95)
    cc$intensity_range    <- cc$intensity_range    %||% c(0.1, 0.5)
    cc$exclusion_eps_rt   <- cc$exclusion_eps_rt   %||% 0.02
    cc$exclusion_eps_cv   <- cc$exclusion_eps_cv   %||% 0.05
    n_con <- if (length(cc$n_per_sample) == 1L) as.integer(cc$n_per_sample)
             else                                as.integer(stats::runif(1,
                                                    cc$n_per_sample[1],
                                                    cc$n_per_sample[2] + 1))
    if (n_con > 0L) {
      srp <- unlist(cmp$sigma_rt_obs); scp <- unlist(cmp$sigma_cv_obs)
      pos_rt <- if (!is.null(cmp$rt_frac)) cmp$rt_frac else cmp$rt_loc
      pos_cv <- if (!is.null(cmp$cv_frac)) cmp$cv_frac else cmp$cv_loc
      keep_pos <- is.finite(pos_rt) & is.finite(pos_cv)
      ocp <- list(
        rt_range = cc$rt_range, cv_range = cc$cv_range,
        intensity_range = cc$intensity_range,
        exclusion_eps_rt = cc$exclusion_eps_rt,
        exclusion_eps_cv = cc$exclusion_eps_cv,
        sigma_rt_pool = srp[is.finite(srp) & srp > 0],
        sigma_cv_pool = scp[is.finite(scp) & scp > 0],
        cat_rt_frac   = pos_rt[keep_pos],
        cat_cv_frac   = pos_cv[keep_pos]
      )
      if (length(ocp$sigma_rt_pool) == 0L ||
          length(ocp$sigma_cv_pool) == 0L) {
        # No usable sigma observations — skip contamination for this
        # sample rather than blowing up the training loop.
        n_con <- 0L
      }
      inj <- inject_off_catalog(Z, ocp, H, W, n_con,
                                  size_jitter,
                                  location_jitter_rt,
                                  location_jitter_cv)
      Z <- inj$Z
      contamination_positions <- inj$positions
    }
  }

  noisy <- Z
  if (add_noise && !is.null(catalog$noise) &&
      is.finite(catalog$noise$mean) && is.finite(catalog$noise$sd)) {
    noise_sd   <- max(1e-8, catalog$noise$sd * noise_scale)
    noise_mean <- catalog$noise$mean
    if (!is.finite(noise_sd))   noise_sd   <- 1e-8
    if (!is.finite(noise_mean)) noise_mean <- 0
    noise <- matrix(abs(stats::rnorm(H * W,
                                      mean = noise_mean,
                                      sd = noise_sd)),
                     nrow = H, ncol = W)
    noisy <- Z + noise
  }
  stopifnot(
    "generate_noisy_pretrain_sample produced non-finite pixels" =
      all(is.finite(Z)) && all(is.finite(noisy)))

  list(clean = Z, noisy = noisy,
       anchor_ids = anchor_ids,
       variable_ids = variable_ids,
       contamination_positions = contamination_positions)
}
