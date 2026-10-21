# classydms: peak catalog construction and diagnostics
# Groups pool peaks into "compounds" — specific (RT, CV) positions
# where peaks reliably appear across samples — for realistic synthetic
# data generation. See PIPELINE_KNOBS.md and the case-control tutorial
# for context.

# ---------------------------------------------------------------------------
# peak_params field registry
#
# Every top-level field ending in `_raw` is a per-peak vector, parallel to
# every other `_raw` field (same length, same peak order). Per-sample fields
# are listed explicitly below. Transforms that drop or reorder peaks /
# samples must go through .subset_peaks() / .subset_samples() so new fields
# that follow these conventions are carried along automatically.
# ---------------------------------------------------------------------------
.pp_per_sample_vectors   <- c("sample_names", "sample_cohort")
.pp_per_sample_matrices  <- c("time_ranges", "cv_ranges", "data_dims")
.pp_alignment_per_sample <- c("shifts_rt", "shifts_cv", "n_anchors_matched")

pp_per_peak_fields <- function(peak_params) {
  nm <- grep("_raw$", names(peak_params), value = TRUE)
  nm[!vapply(peak_params[nm], is.null, logical(1))]
}

#' Validate the internal consistency of a peak_params object
#'
#' Checks the invariants every function that consumes `peak_params`
#' relies on:
#' \itemize{
#'   \item all per-peak vectors (every field ending in `_raw`, e.g.
#'     `rt_loc_raw`, `rt_seconds_raw`, `sample_idx_raw`) have one
#'     common length;
#'   \item `rt_seconds_raw` and `cv_volts_raw` are present together or
#'     not at all;
#'   \item `sample_idx_raw`, if present, is integer-valued with no `NA`
#'     and lies in `1..n_samples`;
#'   \item per-sample fields (`sample_names`, `sample_cohort`,
#'     `n_peaks$values`, the per-sample vectors in `alignment`) have
#'     length `n_samples`, and `time_ranges` / `cv_ranges` have
#'     `n_samples` rows;
#'   \item `rt_loc$values` / `cv_loc$values`, if present, match the
#'     per-peak length.
#' }
#' Called at the top of [build_peak_catalog()] and at the end of every
#' function that returns a transformed `peak_params`. Useful on saved
#' objects too: a physical-mode `peak_params` written by
#' [subset_peak_params()] before this check existed fails here (its
#' `rt_seconds_raw` / `cv_volts_raw` were left unsubset).
#'
#' @param peak_params A `peak_params` list.
#' @param context Label used in error messages.
#' @return `peak_params`, invisibly. Errors on the first violation.
#' @export
validate_peak_params <- function(peak_params, context = "peak_params") {
  fail <- function(...) stop(context, ": ", ..., call. = FALSE)
  if (!is.list(peak_params)) fail("not a list.")
  raw <- pp_per_peak_fields(peak_params)
  if (!length(raw)) fail("has no per-peak (*_raw) vectors.")
  lens <- vapply(peak_params[raw], length, integer(1))
  if (length(unique(lens)) > 1L) {
    fail("per-peak vectors have different lengths (",
         paste(sprintf("%s = %d", raw, lens), collapse = ", "),
         "). They must stay parallel; was this object produced by a ",
         "subset/merge that skipped some fields?")
  }
  n_peaks <- lens[[1]]
  if (xor(is.null(peak_params$rt_seconds_raw),
          is.null(peak_params$cv_volts_raw))) {
    fail("has only one of rt_seconds_raw / cv_volts_raw.")
  }
  for (nm in c("rt_loc", "cv_loc")) {
    v <- peak_params[[nm]]$values
    if (!is.null(v) && length(v) != n_peaks)
      fail(nm, "$values has length ", length(v), ", expected ", n_peaks, ".")
  }

  n_s <- peak_params$n_samples
  if (is.null(n_s)) return(invisible(peak_params))
  if (length(n_s) != 1L || !is.finite(n_s) || n_s < 0 || n_s != round(n_s))
    fail("n_samples must be a single non-negative integer.")
  sidx <- peak_params$sample_idx_raw
  if (!is.null(sidx)) {
    if (anyNA(sidx) || any(sidx != round(sidx)) ||
        any(sidx < 1) || any(sidx > n_s))
      fail("sample_idx_raw must be integers in 1..n_samples (", n_s, ").")
  }
  check_len <- function(x, label) {
    if (!is.null(x) && length(x) != n_s)
      fail(label, " has length ", length(x), ", expected n_samples = ", n_s, ".")
  }
  for (nm in .pp_per_sample_vectors) check_len(peak_params[[nm]], nm)
  check_len(peak_params[["n_peaks"]][["values"]], "n_peaks$values")
  for (nm in .pp_alignment_per_sample)
    check_len(peak_params$alignment[[nm]], paste0("alignment$", nm))
  for (nm in .pp_per_sample_matrices) {
    m <- peak_params[[nm]]
    if (!is.null(m) && NROW(m) != n_s)
      fail(nm, " has ", NROW(m), " rows, expected n_samples = ", n_s, ".")
  }
  invisible(peak_params)
}

# Keep peaks `keep` (logical or integer index) in every per-peak vector and
# refresh the fields derived from them.
.subset_peaks <- function(peak_params, keep) {
  for (nm in pp_per_peak_fields(peak_params))
    peak_params[[nm]] <- peak_params[[nm]][keep]
  .refresh_location_summaries(peak_params)
}

.refresh_location_summaries <- function(peak_params) {
  loc <- function(x) list(
    values = x,
    mean   = if (length(x)) mean(x, na.rm = TRUE) else NA_real_,
    sd     = if (length(x) > 1L) stats::sd(x, na.rm = TRUE) else 0)
  if (!is.null(peak_params$rt_loc_raw)) peak_params$rt_loc <- loc(peak_params$rt_loc_raw)
  if (!is.null(peak_params$cv_loc_raw)) peak_params$cv_loc <- loc(peak_params$cv_loc_raw)
  peak_params$n_peaks_detected <- length(peak_params$rt_loc_raw)
  peak_params
}

# Keep samples `idx` (integer, in order) in every per-sample field.
.subset_samples <- function(peak_params, idx) {
  for (nm in .pp_per_sample_vectors)
    if (!is.null(peak_params[[nm]])) peak_params[[nm]] <- peak_params[[nm]][idx]
  for (nm in .pp_per_sample_matrices)
    if (!is.null(peak_params[[nm]]))
      peak_params[[nm]] <- peak_params[[nm]][idx, , drop = FALSE]
  if (!is.null(peak_params$alignment)) {
    for (nm in .pp_alignment_per_sample)
      if (!is.null(peak_params$alignment[[nm]]))
        peak_params$alignment[[nm]] <- peak_params$alignment[[nm]][idx]
  }
  if (!is.null(peak_params[["n_peaks"]][["values"]])) {
    v <- peak_params[["n_peaks"]][["values"]][idx]
    peak_params$n_peaks$values <- v
    peak_params$n_peaks$mean <- if (length(v)) mean(v) else NA_real_
    peak_params$n_peaks$sd   <- if (length(v) > 1L) stats::sd(v) else 0
  }
  peak_params$n_samples <- length(idx)
  peak_params
}

.refit_lognormal_summaries <- function(peak_params) {
  fit_lognormal <- function(x) {
    x <- x[is.finite(x) & x > 0]
    if (length(x) < 10L) return(list(meanlog = 0, sdlog = 1))
    list(meanlog = mean(log(x)), sdlog = stats::sd(log(x)))
  }
  if (!is.null(peak_params$sigma_rt_raw)) {
    peak_params$sigma_rt  <- fit_lognormal(peak_params$sigma_rt_raw)
    peak_params$sigma_cv  <- fit_lognormal(peak_params$sigma_cv_raw)
    peak_params$intensity <- fit_lognormal(peak_params$intensity_raw)
  }
  peak_params
}

# Default loose eps for the alignment pre-pass. Alignment always clusters in
# fraction space, so in physical mode the catalog eps (seconds / volts) is
# converted with the median per-sample axis span.
.default_align_eps <- function(eps, axis_ranges, coord_mode, fn) {
  if (coord_mode == "fraction") return(3 * eps)
  if (is.null(axis_ranges))
    stop(fn, ": align = TRUE with coord_mode = 'physical' needs ",
         "time_ranges / cv_ranges on peak_params to convert eps to ",
         "fractions; pass eps_rt_loose / eps_cv_loose explicitly.",
         call. = FALSE)
  3 * eps / stats::median(axis_ranges[, "max"] - axis_ranges[, "min"])
}

#' Backfill sample_idx_raw and sample_names on an older peak_params
#'
#' `sample_idx_raw` and `sample_names` were added to
#' [estimate_peak_params()] in a recent commit. `peak_params` objects
#' saved before then lack these fields and cannot be used with
#' [build_peak_catalog()] directly. This helper reconstructs
#' `sample_idx_raw` from the per-sample peak count already stored in
#' `n_peaks$values` and optionally attaches `sample_names` from a
#' user-supplied vector.
#'
#' Peaks in the `_raw` vectors are stored in sample-then-detection
#' order (all of sample 1's peaks, then all of sample 2's peaks, ...),
#' and each sample contributes exactly `n_peaks$values[i]` entries.
#' The reconstruction is `rep(1:n_samples, n_peaks$values)`.
#'
#' Sanity-checks that `sum(n_peaks$values) == length(rt_loc_raw)`; if
#' this fails, the peak_params object is inconsistent and must be
#' regenerated by rerunning `estimate_peak_params`.
#'
#' @param peak_params An older `peak_params` object missing
#'   `sample_idx_raw`.
#' @param sample_names Optional character vector of sample names.
#'   Must have length `n_samples`. If `NULL`, `sample_names` is not
#'   attached (the catalog doesn't need it; only the peak-params
#'   report does).
#' @return A `peak_params`-shaped list with `sample_idx_raw` and
#'   (if provided) `sample_names` added.
#' @export
backfill_sample_idx_raw <- function(peak_params, sample_names = NULL) {
  if (is.null(peak_params[["n_peaks"]][["values"]])) {
    stop("backfill_sample_idx_raw: peak_params is missing ",
         "n_peaks$values — cannot reconstruct sample indices. ",
         "Rerun estimate_peak_params from Z data.")
  }
  if (is.null(peak_params$rt_loc_raw)) {
    stop("backfill_sample_idx_raw: peak_params is missing rt_loc_raw. ",
         "This object predates the empirical-location commit. ",
         "Rerun estimate_peak_params.")
  }
  n_peaks_per_sample <- peak_params[["n_peaks"]][["values"]]
  expected_total <- sum(n_peaks_per_sample)
  observed_total <- length(peak_params$rt_loc_raw)
  if (expected_total != observed_total) {
    stop("backfill_sample_idx_raw: inconsistent peak counts — ",
         "sum(n_peaks$values) = ", expected_total,
         " but length(rt_loc_raw) = ", observed_total,
         ". Rerun estimate_peak_params from Z data.")
  }
  peak_params$sample_idx_raw <-
    rep(seq_along(n_peaks_per_sample), times = n_peaks_per_sample)
  if (is.null(peak_params$n_samples)) {
    peak_params$n_samples <- length(n_peaks_per_sample)
  }
  if (!is.null(sample_names)) {
    if (length(sample_names) != peak_params$n_samples) {
      stop("backfill_sample_idx_raw: sample_names must have length ",
           peak_params$n_samples, ", got ", length(sample_names))
    }
    peak_params$sample_names <- sample_names
  }
  validate_peak_params(peak_params, "backfill_sample_idx_raw")
}

#' Subset a peak_params object to a chosen set of samples
#'
#' Given a `peak_params` object with `sample_idx_raw` (or after
#' [backfill_sample_idx_raw()]), returns a new `peak_params` containing
#' only the peaks and metadata for the specified samples. Useful for
#' building per-condition catalogs from a single pooled
#' `estimate_peak_params` run — split a combined disease+control
#' `peak_params` into `pp_disease` and `pp_control` without rerunning
#' the expensive step.
#'
#' Sample indices in the returned object are remapped to `1:length(sample_indices)`
#' so downstream tools that expect contiguous 1-based indexing work.
#' Samples keep their original relative order (`sample_indices` is
#' sorted). Every per-peak vector (every field ending in `_raw`,
#' including `rt_seconds_raw` / `cv_volts_raw`) is subset with the same
#' mask, and every per-sample field (`sample_names`, `sample_cohort`,
#' `time_ranges`, `cv_ranges`, `n_peaks$values`, per-sample `alignment`
#' vectors) is subset to the kept samples. Fitted summary
#' statistics (`rt_loc$mean/sd`, `sigma_rt$meanlog/sdlog`, etc.) are
#' recomputed from the subset; noise stats are kept as-is (pooled
#' across samples in the original run, which is usually close enough
#' unless the condition of interest has systematically different
#' background noise).
#'
#' @param peak_params A `peak_params` object with `sample_idx_raw`.
#' @param sample_indices Integer vector of sample indices (in the
#'   original peak_params' 1..n_samples space) to keep. Common
#'   pattern: `which(y == disease_label)`.
#' @return A `peak_params`-shaped list containing only the specified
#'   samples' peaks. Ready to pass to [build_peak_catalog()] and other
#'   downstream tools.
#' @export
subset_peak_params <- function(peak_params, sample_indices) {
  if (is.null(peak_params$sample_idx_raw)) {
    stop("subset_peak_params: peak_params$sample_idx_raw is missing. ",
         "Call backfill_sample_idx_raw() first, or rerun ",
         "estimate_peak_params.")
  }
  validate_peak_params(peak_params, "subset_peak_params (input)")
  sample_indices <- as.integer(sample_indices)
  if (any(is.na(sample_indices)) || any(sample_indices < 1L)) {
    stop("subset_peak_params: sample_indices must be positive integers.")
  }
  if (anyDuplicated(sample_indices)) {
    stop("subset_peak_params: sample_indices contains duplicate values.")
  }
  n_orig <- if (!is.null(peak_params$n_samples)) peak_params$n_samples
            else max(peak_params$sample_idx_raw)
  if (any(sample_indices > n_orig)) {
    stop("subset_peak_params: sample_indices contains values beyond ",
         "the original n_samples (", n_orig, ").")
  }
  kept <- sort(sample_indices)

  out <- .subset_peaks(peak_params, peak_params$sample_idx_raw %in% kept)
  # Remap sample_idx_raw to contiguous 1..length(sample_indices)
  out$sample_idx_raw <- match(out$sample_idx_raw, kept)
  out <- .subset_samples(out, kept)
  out <- .refit_lognormal_summaries(out)

  # noise stays as-is (pooled from original run — usually close enough)

  validate_peak_params(out, "subset_peak_params")
}

#' Build a peak catalog from a peak_params object
#'
#' Greedy clustering: peaks are sorted by intensity descending, each
#' unassigned peak becomes a new cluster seed, and all remaining
#' unassigned peaks within `(eps_rt, eps_cv)` of the seed are added
#' to that cluster. Clusters are then classified into three groups:
#' \describe{
#'   \item{compounds}{clusters with `>= min_cluster_size` observations
#'     AND spanning `>= min_support_frac` of samples. These are the
#'     "reliable positions" the catalog is built around.}
#'   \item{small_clusters}{clusters with >= 2 observations but below
#'     the sample-support threshold. Might be sporadic contamination,
#'     small compounds, or peaks from just a few subjects.}
#'   \item{singletons}{clusters of size 1 — peaks that didn't co-locate
#'     with any other pool peak. Often noise, but worth inspecting.}
#' }
#'
#' Uses a sorted-index optimization so runtime stays `O(n log n + n·k)`
#' rather than `O(n^2)`, where `k` is the average number of candidates
#' per seed. Practical for pool sizes up to millions.
#'
#' @param peak_params Output of [estimate_peak_params()] containing
#'   `rt_loc_raw`, `cv_loc_raw`, `sigma_rt_raw`, `sigma_cv_raw`,
#'   `intensity_raw`, `sample_idx_raw`, and `n_samples`.
#' @param eps_rt,eps_cv Clustering radius in fraction-of-image units.
#'   Defaults `(0.01, 0.05)` roughly match the peak-detector's
#'   `min_sep_rt = 8, min_sep_cv = 2` on a 1400 x 40 native image.
#'   Tighter → more, smaller clusters. Looser → fewer, larger clusters.
#'   In `coord_mode = "physical"` they are in seconds and volts.
#' @param coord_mode `"fraction"` (default) clusters on
#'   `rt_loc_raw` / `cv_loc_raw`; `"physical"` clusters on
#'   `rt_seconds_raw` / `cv_volts_raw` (requires [estimate_peak_params()]
#'   to have been given `time_axes` / `cv_axes`). Either way the catalog
#'   stores fraction positions (`rt_frac`, `cv_frac`, `*_frac_obs`) and,
#'   when available, physical ones (`rt_seconds`, `cv_volts`,
#'   `rt_seconds_obs`, `cv_volts_obs`), all from the same member peaks.
#' @param min_support_frac Fraction of samples a cluster must span to
#'   count as a compound. Default `0.10`: appear in at least 10% of
#'   samples.
#' @param min_cluster_size Minimum number of pool peaks in a cluster
#'   for compound status. Default `2L` (excludes singletons from
#'   compounds regardless of support).
#'
#' @return An object of class `peak_catalog` with fields:
#' \describe{
#'   \item{compounds}{tibble of compound entries — one row per
#'     compound with columns `compound_id, rt_loc, cv_loc,
#'     location_sd_rt, location_sd_cv, n_observations, n_samples,
#'     prevalence`, and list-columns `sigma_rt_obs, sigma_cv_obs,
#'     intensity_obs, rt_frac_obs, cv_frac_obs, member_indices`, plus
#'     `rt_seconds_obs, cv_volts_obs` when physical positions exist and
#'     `sample_idx_obs` (source sample of each observation) when data
#'     dims are known.}
#'   \item{small_clusters}{same schema, for clusters below thresholds.}
#'   \item{singletons}{same schema, for clusters of size 1.}
#'   \item{parameters}{the settings used to build the catalog, plus
#'     `data_dims`: the `n_samples x 2` matrix of each sample's unpadded
#'     `rows, cols` (from `peak_params$data_dims`, or the axis lengths
#'     in `time_ranges` / `cv_ranges`), which generators use to place
#'     observations on a centre-padded canvas
#'     ([catalog_pixel_positions()]).}
#' }
#' @param align Optional per-sample drift correction applied to
#'   `peak_params` before clustering. `FALSE` (default) skips
#'   alignment. `TRUE` calls [align_via_catalog()] with sensible
#'   defaults (`eps_rt_loose = 3 * eps_rt`, `eps_cv_loose = 3 * eps_cv`,
#'   `min_support_frac = min_support_frac`, `iterate = TRUE`,
#'   `max_iter = 3L`). Pass a named list to override any of those.
#'   Alignment is idempotent and safe to leave off if drift is small;
#'   see the compactness-delta table in
#'   `compare_all_catalogs.Rmd` for how to decide whether it helps
#'   on your cohort.
#' @export
build_peak_catalog <- function(peak_params,
                                eps_rt = 0.01, eps_cv = 0.05,
                                min_support_frac = 0.10,
                                min_cluster_size = 2L,
                                coord_mode = c("fraction", "physical"),
                                align = FALSE) {
  coord_mode <- match.arg(coord_mode)
  validate_peak_params(peak_params, "build_peak_catalog")
  if (!isFALSE(align)) {
    align_args <- if (isTRUE(align)) list() else align
    if (!is.list(align_args))
      stop("build_peak_catalog: `align` must be FALSE, TRUE, or a ",
           "named list of align_via_catalog parameters.")
    align_args$eps_rt_loose     <- align_args$eps_rt_loose     %||%
      .default_align_eps(eps_rt, peak_params$time_ranges, coord_mode,
                         "build_peak_catalog")
    align_args$eps_cv_loose     <- align_args$eps_cv_loose     %||%
      .default_align_eps(eps_cv, peak_params$cv_ranges, coord_mode,
                         "build_peak_catalog")
    align_args$min_support_frac <- align_args$min_support_frac %||% min_support_frac
    align_args$iterate          <- align_args$iterate          %||% TRUE
    align_args$max_iter         <- align_args$max_iter         %||% 3L
    align_args$verbose          <- align_args$verbose          %||% FALSE
    message("build_peak_catalog: running align_via_catalog ",
            "(eps_rt_loose = ", align_args$eps_rt_loose,
            ", eps_cv_loose = ", align_args$eps_cv_loose, ") ...")
    peak_params <- do.call(align_via_catalog,
                            c(list(peak_params = peak_params), align_args))
  }
  needed <- c("rt_loc_raw", "cv_loc_raw", "sample_idx_raw",
              "sigma_rt_raw", "sigma_cv_raw", "intensity_raw",
              "n_samples")
  missing <- needed[vapply(needed, function(nm) is.null(peak_params[[nm]]),
                            logical(1))]
  if (length(missing)) {
    stop("build_peak_catalog: peak_params is missing ",
         paste(missing, collapse = ", "), ". Ensure sample_idx_raw ",
         "was added by a recent version of estimate_peak_params().")
  }
  if (coord_mode == "physical") {
    if (is.null(peak_params$rt_seconds_raw) ||
        is.null(peak_params$cv_volts_raw)) {
      stop("coord_mode = 'physical' requires peak_params to have ",
           "rt_seconds_raw and cv_volts_raw fields. Re-run ",
           "estimate_peak_params() passing time_axes and cv_axes.")
    }
  }

  # Choose coordinate vectors: physical (seconds / volts) or fraction.
  # Only the clustering itself is affected — the returned catalog stores
  # BOTH sets so downstream tools can use whichever they want.
  rt <- if (coord_mode == "physical") peak_params$rt_seconds_raw
        else                          peak_params$rt_loc_raw
  cv <- if (coord_mode == "physical") peak_params$cv_volts_raw
        else                          peak_params$cv_loc_raw
  # Preserve the alternate representation for storage
  rt_frac_ref <- peak_params$rt_loc_raw
  cv_frac_ref <- peak_params$cv_loc_raw
  intensity  <- peak_params$intensity_raw
  sample_idx <- peak_params$sample_idx_raw
  n_peaks    <- length(rt)
  n_samples  <- peak_params$n_samples

  if (n_peaks == 0L) {
    stop("build_peak_catalog: peak_params$rt_loc_raw is empty.")
  }
  # Every vector indexed by cluster members must share the clustering
  # vector's length; validate_peak_params() guarantees this, but the
  # members are applied to all of them, so assert it where it matters.
  stopifnot(length(cv) == n_peaks, length(rt_frac_ref) == n_peaks,
            length(cv_frac_ref) == n_peaks, length(intensity) == n_peaks,
            length(sample_idx) == n_peaks)
  rt_phys <- peak_params$rt_seconds_raw
  cv_phys <- peak_params$cv_volts_raw
  # Unpadded per-sample data dims, so generators can place each
  # observation on a centre-padded canvas (see R/placement.R).
  data_dims <- .pp_data_dims(peak_params)

  # Sort peaks by rt_loc for fast range lookup (findInterval).
  rt_order  <- order(rt)
  rt_sorted <- rt[rt_order]

  # Greedy assignment in intensity-descending order.
  intensity_order <- order(intensity, decreasing = TRUE)
  cluster_id <- rep(NA_integer_, n_peaks)
  current_cluster <- 0L

  for (i in intensity_order) {
    if (!is.na(cluster_id[i])) next
    current_cluster <- current_cluster + 1L
    cluster_id[i]   <- current_cluster

    # Binary search for peaks whose rt is within [rt[i] - eps_rt, rt[i] + eps_rt].
    lo <- findInterval(rt[i] - eps_rt, rt_sorted, all.inside = FALSE) + 1L
    hi <- findInterval(rt[i] + eps_rt, rt_sorted, all.inside = FALSE)
    if (hi < lo) next
    candidates <- rt_order[lo:hi]
    # Filter by cv proximity and unassigned status.
    keep <- is.na(cluster_id[candidates]) &
            abs(cv[candidates] - cv[i]) <= eps_cv
    if (!any(keep)) next
    cluster_id[candidates[keep]] <- current_cluster
  }

  # Build a per-cluster summary. Use vapply for speed and tibble
  # with list-columns for downstream ergonomics.
  n_clusters <- current_cluster

  # Precompute per-cluster member index vectors.
  cluster_members <- split(seq_len(n_peaks), cluster_id)

  n_cohort_samples <- n_samples   # avoid shadowing inside tibble()
  cluster_rows <- lapply(seq_len(n_clusters), function(cid) {
    members <- cluster_members[[as.character(cid)]]
    unique_samples <- unique(sample_idx[members])
    n_uniq <- length(unique_samples)
    # Build the base row using the CLUSTERING coordinate system as
    # rt_loc / cv_loc (fraction OR physical depending on coord_mode).
    # NA-safe aggregation across cluster members. Raw peak_params
    # occasionally contain NAs on individual axes (e.g., degenerate
    # peaks near the pad boundary); without na.rm the whole compound's
    # rt_loc / rt_frac / rt_seconds becomes NA and every downstream
    # generation of that compound produces NaN pixels. Also filter the
    # per-observation morphology list-columns so downstream bootstrap
    # draws can't pick an NA.
    sig_rt_m <- peak_params$sigma_rt_raw[members]
    sig_cv_m <- peak_params$sigma_cv_raw[members]
    int_m    <- intensity[members]
    # Per-observation fraction positions of each cluster member.
    # Stored as list-columns so synthesis can bootstrap position jointly
    # with morphology (same random index j across all four vectors), which
    # gives synthetic peaks a natural spread ~= location_sd_rt/cv rather
    # than the tiny location_jitter that was previously the only source
    # of positional variability. Preserves the physical coupling that
    # peaks with certain morphologies also land at certain positions.
    rt_frac_m <- rt_frac_ref[members]
    cv_frac_m <- cv_frac_ref[members]
    keep_obs <- is.finite(sig_rt_m) & sig_rt_m > 0 &
                is.finite(sig_cv_m) & sig_cv_m > 0 &
                is.finite(int_m)    & int_m    > 0 &
                is.finite(rt_frac_m) & is.finite(cv_frac_m)
    row <- tibble::tibble(
      compound_id     = cid,
      rt_loc          = stats::median(rt[members], na.rm = TRUE),
      cv_loc          = stats::median(cv[members], na.rm = TRUE),
      location_sd_rt  = if (length(members) > 1L)
                          stats::sd(rt[members], na.rm = TRUE) else 0,
      location_sd_cv  = if (length(members) > 1L)
                          stats::sd(cv[members], na.rm = TRUE) else 0,
      n_observations  = length(members),
      n_samples       = n_uniq,
      prevalence      = n_uniq / n_cohort_samples,
      sigma_rt_obs    = list(sig_rt_m[keep_obs]),
      sigma_cv_obs    = list(sig_cv_m[keep_obs]),
      intensity_obs   = list(int_m[keep_obs]),
      rt_frac_obs     = list(rt_frac_m[keep_obs]),
      cv_frac_obs     = list(cv_frac_m[keep_obs]),
      member_indices  = list(members)
    )
    # Per-observation physical positions (same keep_obs filter), so
    # downstream code can map observations onto a target pixel grid via
    # the instrument axes instead of per-sample fractions.
    if (!is.null(data_dims))
      row$sample_idx_obs <- list(as.integer(sample_idx[members][keep_obs]))
    if (!is.null(rt_phys)) {
      row$rt_seconds_obs <- list(rt_phys[members][keep_obs])
      row$cv_volts_obs   <- list(cv_phys[members][keep_obs])
    }
    # If both physical + fraction are available, also store the OTHER
    # representation as rt_seconds / cv_volts (or rt_frac / cv_frac).
    if (coord_mode == "physical") {
      row$rt_seconds <- row$rt_loc   # already in seconds
      row$cv_volts   <- row$cv_loc   # already in volts
      row$rt_frac    <- stats::median(rt_frac_ref[members], na.rm = TRUE)
      row$cv_frac    <- stats::median(cv_frac_ref[members], na.rm = TRUE)
    } else if (!is.null(peak_params$rt_seconds_raw)) {
      row$rt_frac    <- row$rt_loc   # already in fraction
      row$cv_frac    <- row$cv_loc   # already in fraction
      row$rt_seconds <- stats::median(rt_phys[members], na.rm = TRUE)
      row$cv_volts   <- stats::median(cv_phys[members], na.rm = TRUE)
    }
    row
  })
  clusters_df <- dplyr::bind_rows(cluster_rows)

  # Members whose morphology / fraction position is unusable are left out
  # of the *_obs list-columns (but still count toward n_observations and
  # prevalence). Report how many, rather than dropping them silently.
  n_obs_invalid <- sum(clusters_df$n_observations) -
                   sum(lengths(clusters_df$intensity_obs))
  if (n_obs_invalid > 0L) {
    message("build_peak_catalog: ", n_obs_invalid, " of ", n_peaks,
            " pool peaks have non-finite or non-positive sigma / ",
            "intensity / fraction position and are excluded from the ",
            "*_obs list-columns.")
  }

  # Classify into compounds / small_clusters / singletons.
  # Reject any cluster with a non-finite coordinate OR no usable
  # per-observation morphology draws (empty intensity_obs list) — the
  # generator can't render these safely.
  n_obs_valid <- vapply(clusters_df$intensity_obs, length, integer(1))
  is_bad_coord <- !is.finite(clusters_df$rt_loc) |
                  !is.finite(clusters_df$cv_loc) |
                  n_obs_valid == 0L
  if ("rt_frac" %in% names(clusters_df))
    is_bad_coord <- is_bad_coord | !is.finite(clusters_df$rt_frac) |
                                    !is.finite(clusters_df$cv_frac)
  if (any(is_bad_coord)) {
    message("build_peak_catalog: dropping ", sum(is_bad_coord),
            " cluster(s) with non-finite coordinates or no usable ",
            "per-observation morphology.")
  }
  is_singleton      <- clusters_df$n_observations < 2L
  is_below_support  <- clusters_df$prevalence < min_support_frac
  is_below_size     <- clusters_df$n_observations < min_cluster_size
  is_compound       <- !is_singleton & !is_below_support & !is_below_size &
                       !is_bad_coord

  compounds       <- clusters_df[is_compound, ]
  small_clusters  <- clusters_df[!is_compound & !is_singleton & !is_bad_coord, ]
  singletons      <- clusters_df[is_singleton & !is_bad_coord, ]

  # Re-number compounds starting from 1 for convenience.
  if (nrow(compounds) > 0L) {
    compounds$compound_id <- seq_len(nrow(compounds))
  }

  structure(
    list(
      compounds       = compounds,
      small_clusters  = small_clusters,
      singletons      = singletons,
      noise           = if (!is.null(peak_params$noise))
                          peak_params$noise
                        else list(mean = 0.005, sd = 0.005),
      parameters      = list(
        eps_rt              = eps_rt,
        eps_cv              = eps_cv,
        min_support_frac    = min_support_frac,
        min_cluster_size    = min_cluster_size,
        n_pool_peaks        = n_peaks,
        n_obs_invalid       = n_obs_invalid,
        n_samples           = n_samples,
        n_total_clusters    = n_clusters,
        coord_mode          = coord_mode,
        data_dims           = data_dims
      )
    ),
    class = "peak_catalog"
  )
}

#' Diagnostic summary of a peak catalog
#'
#' Returns a list of scalar and small-data-frame summary stats and,
#' when `print = TRUE`, prints a human-readable summary. Intended as
#' the first thing you look at after `build_peak_catalog()`.
#'
#' @param catalog A `peak_catalog` object.
#' @param prevalence_thresholds Numeric vector of prevalence values at
#'   which to report the number of compounds meeting or exceeding that
#'   threshold. Default `c(0.10, 0.25, 0.50, 0.75, 0.90, 0.99)`.
#' @param print If `TRUE` (default), print a formatted summary in
#'   addition to returning the underlying data.
#'
#' @return List with `n_pool_peaks`, `n_samples`, `n_compounds`,
#'   `n_small_clusters`, `n_singletons`, `total_observations`,
#'   `fraction_singleton`, `fraction_below_support`, and a data
#'   frame `compounds_at_threshold` with columns
#'   `prevalence_threshold, n_compounds`.
#' @export
catalog_summary <- function(catalog,
                             prevalence_thresholds = c(0.10, 0.25, 0.50,
                                                         0.75, 0.90, 0.99),
                             print = TRUE) {
  stopifnot(inherits(catalog, "peak_catalog"))

  n_pool   <- catalog$parameters$n_pool_peaks
  n_samp   <- catalog$parameters$n_samples
  n_comp   <- nrow(catalog$compounds)
  n_small  <- nrow(catalog$small_clusters)
  n_singl  <- nrow(catalog$singletons)
  n_total  <- n_comp + n_small + n_singl

  # Total observations across all cluster types
  count_obs <- function(df) if (nrow(df) == 0L) 0L
                            else sum(df$n_observations)
  total_obs <- count_obs(catalog$compounds) +
               count_obs(catalog$small_clusters) +
               count_obs(catalog$singletons)

  # Prevalence sweep — but the below-support clusters have low prevalence
  # by definition, so this is over compounds only.
  n_at <- vapply(prevalence_thresholds, function(p)
    sum(catalog$compounds$prevalence >= p), integer(1))
  compounds_at_threshold <- tibble::tibble(
    prevalence_threshold = prevalence_thresholds,
    n_compounds          = n_at
  )

  out <- list(
    n_pool_peaks            = n_pool,
    n_samples               = n_samp,
    n_total_clusters        = n_total,
    n_compounds             = n_comp,
    n_small_clusters        = n_small,
    n_singletons            = n_singl,
    total_observations      = total_obs,
    fraction_singleton      = if (n_total) n_singl / n_total else NA_real_,
    fraction_below_support  = if (n_total) n_small / n_total else NA_real_,
    compounds_at_threshold  = compounds_at_threshold
  )

  if (isTRUE(print)) {
    cat(sprintf("Peak catalog summary\n"))
    cat(sprintf("--------------------\n"))
    cat(sprintf("  Pool peaks:            %s\n",
                format(n_pool, big.mark = ",")))
    cat(sprintf("  Cohort samples:        %d\n", n_samp))
    cat(sprintf("  Clusters found:        %s\n",
                format(n_total, big.mark = ",")))
    cat(sprintf("    compounds (support >= %.0f%% and size >= %d): %d\n",
                100 * catalog$parameters$min_support_frac,
                catalog$parameters$min_cluster_size, n_comp))
    cat(sprintf("    small_clusters (size >= 2 but below support):  %d\n",
                n_small))
    cat(sprintf("    singletons (size == 1):                        %d\n",
                n_singl))
    cat(sprintf("  Fraction singletons:   %.1f%%\n",
                100 * out$fraction_singleton))
    cat(sprintf("  Fraction below support:%.1f%%\n",
                100 * out$fraction_below_support))
    cat(sprintf("\nCompounds at prevalence thresholds:\n"))
    print(compounds_at_threshold, n = Inf)
  }

  invisible(out)
}

#' Compactness diagnostic for a peak catalog
#'
#' Reports statistics that reveal whether the clustering radius is
#' set correctly. The key statistic is `median_obs_per_sample` — the
#' typical number of pool peaks each compound catches per sample. A
#' well-clustered catalog has this near 1.0 (each real compound
#' contributes one peak per sample on average). Values noticeably
#' above 1.0 mean the clusters are wide enough to catch multiple
#' distinct peaks from the same sample, indicating over-merging — a
#' looser eps than the true compound-position separation.
#'
#' Also reports the distribution of location-spread within compounds
#' (in pixel units at a specified target image size, default native
#' `1400 x 40`), and the fraction of compounds passing progressively
#' stricter prevalence thresholds.
#'
#' @param catalog A `peak_catalog` object.
#' @param H,W Image dimensions in pixels used to convert location
#'   spread from fraction back to pixels. Defaults `1400 x 40` (native
#'   GC-DMS). Pass your actual padded / trimmed dimensions for
#'   accurate diagnostics.
#' @param over_merge_threshold Cutoff above which a compound is
#'   flagged as over-merged. Default `1.5` (a compound catching more
#'   than 1.5 pool peaks per sample on average is suspect).
#'
#' @return A one-row tibble with compactness statistics. Also invisibly
#'   returns the underlying per-compound data if you want to inspect
#'   individual compounds via `attr(x, "per_compound")`.
#' @export
catalog_compactness <- function(catalog,
                                 H = 1400L, W = 40L,
                                 over_merge_threshold = 1.5) {
  stopifnot(inherits(catalog, "peak_catalog"))
  if (nrow(catalog$compounds) == 0L) {
    return(tibble::tibble(
      n_compounds            = 0L,
      median_obs_per_sample  = NA_real_,
      q75_obs_per_sample     = NA_real_,
      max_obs_per_sample     = NA_real_,
      frac_over_merged       = NA_real_,
      median_loc_sd_rt_px    = NA_real_,
      q95_loc_sd_rt_px       = NA_real_,
      median_loc_sd_cv_px    = NA_real_,
      q95_loc_sd_cv_px       = NA_real_,
      pct_anchor_gte_90pct   = NA_real_,
      pct_anchor_gte_99pct   = NA_real_
    ))
  }
  per_compound <- dplyr::mutate(catalog$compounds,
    obs_per_sample = n_observations / n_samples,
    loc_sd_rt_px   = location_sd_rt * H,
    loc_sd_cv_px   = location_sd_cv * W
  )
  qtile <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE))
  out <- tibble::tibble(
    n_compounds            = nrow(catalog$compounds),
    median_obs_per_sample  = stats::median(per_compound$obs_per_sample),
    q75_obs_per_sample     = qtile(per_compound$obs_per_sample, 0.75),
    max_obs_per_sample     = max(per_compound$obs_per_sample),
    frac_over_merged       = mean(per_compound$obs_per_sample >
                                     over_merge_threshold),
    median_loc_sd_rt_px    = stats::median(per_compound$loc_sd_rt_px),
    q95_loc_sd_rt_px       = qtile(per_compound$loc_sd_rt_px, 0.95),
    median_loc_sd_cv_px    = stats::median(per_compound$loc_sd_cv_px),
    q95_loc_sd_cv_px       = qtile(per_compound$loc_sd_cv_px, 0.95),
    pct_anchor_gte_90pct   = 100 * mean(catalog$compounds$prevalence >= 0.90),
    pct_anchor_gte_99pct   = 100 * mean(catalog$compounds$prevalence >= 0.99)
  )
  attr(out, "per_compound") <- per_compound
  out
}

#' Merge multiple peak_params objects into a single super-cohort
#'
#' Concatenates all `_raw` vectors across the input peak_params
#' objects with globally-unique `sample_idx_raw`. Useful for building
#' a universal catalog from all cohorts jointly.
#'
#' Each returned merged sample's original cohort is recorded in the
#' `sample_cohort` field, which is a character vector of length
#' `n_samples` giving the input list's name for each sample.
#'
#' @param pp_list Named list of `peak_params` objects. All must have
#'   `sample_idx_raw`. Names become cohort labels.
#' @return A single `peak_params`-shaped list with concatenated
#'   `_raw` vectors, renumbered `sample_idx_raw`, `sample_cohort`
#'   giving per-sample cohort membership, and `cohort_names` (input
#'   names).
#' @export
merge_peak_params <- function(pp_list) {
  stopifnot(is.list(pp_list), length(pp_list) >= 1L)
  if (is.null(names(pp_list)))
    names(pp_list) <- paste0("cohort_", seq_along(pp_list))
  for (i in seq_along(pp_list)) {
    if (is.null(pp_list[[i]]$sample_idx_raw)) {
      stop("pp_list[[", i, "]] (", names(pp_list)[i],
           ") is missing sample_idx_raw. ",
           "Call backfill_sample_idx_raw() first.")
    }
  }

  for (i in seq_along(pp_list))
    validate_peak_params(pp_list[[i]],
                         paste0("merge_peak_params: pp_list$", names(pp_list)[i]))

  cohort_n_samples <- vapply(pp_list, function(pp) as.integer(pp$n_samples),
                             integer(1))
  cohort_offsets   <- c(0L, cumsum(cohort_n_samples)[-length(pp_list)])
  total_n_samples  <- sum(cohort_n_samples)

  # Concatenate every per-peak field that all inputs carry (including the
  # physical rt_seconds_raw / cv_volts_raw). Fields only some inputs have
  # cannot be kept parallel, so they are dropped with a message.
  field_sets <- lapply(pp_list, pp_per_peak_fields)
  common  <- Reduce(intersect, field_sets)
  dropped <- setdiff(Reduce(union, field_sets), common)
  if (length(dropped)) {
    message("merge_peak_params: dropping per-peak field(s) not present in ",
            "every input: ", paste(dropped, collapse = ", "), ".")
  }
  merged <- list(n_samples = total_n_samples, noise = pp_list[[1]]$noise)
  for (nm in setdiff(common, "sample_idx_raw"))
    merged[[nm]] <- unlist(lapply(pp_list, `[[`, nm), use.names = FALSE)
  # Global sample indices
  merged$sample_idx_raw <- unlist(lapply(seq_along(pp_list), function(i)
    as.integer(pp_list[[i]]$sample_idx_raw) + cohort_offsets[i]))

  # Per-sample fields that every input carries
  all_have <- function(nm) all(vapply(pp_list, function(pp) !is.null(pp[[nm]]),
                                      logical(1)))
  if (all_have("sample_names"))
    merged$sample_names <- unlist(lapply(pp_list, `[[`, "sample_names"),
                                  use.names = FALSE)
  for (nm in .pp_per_sample_matrices)
    if (all_have(nm)) merged[[nm]] <- do.call(rbind, lapply(pp_list, `[[`, nm))
  if (!is.null(merged$rt_seconds_raw)) merged$coord_mode <- "physical"
  merged <- .refresh_location_summaries(merged)
  merged <- .refit_lognormal_summaries(merged)

  # Per-sample counts
  merged$n_peaks <- list(
    values = as.integer(unlist(lapply(pp_list, function(pp) {
      v <- pp[["n_peaks"]][["values"]]
      if (!is.null(v)) v
      else tabulate(pp$sample_idx_raw, nbins = pp$n_samples)
    })))
  )
  merged$n_peaks$mean <- mean(merged$n_peaks$values)
  merged$n_peaks$sd   <- if (length(merged$n_peaks$values) > 1L)
                          stats::sd(merged$n_peaks$values) else 0

  # Cohort membership: which cohort does each sample belong to?
  merged$sample_cohort <- unlist(mapply(
    function(name, n) rep(name, n),
    names(pp_list), cohort_n_samples,
    SIMPLIFY = FALSE), use.names = FALSE)
  merged$cohort_names <- names(pp_list)

  validate_peak_params(merged, "merge_peak_params")
}

#' Build a universal catalog by pooling all samples from multiple cohorts
#'
#' Unlike [merge_catalogs()], which reconciles already-built per-cohort
#' catalogs, this function pools raw peaks from all cohorts BEFORE
#' clustering. Every peak from every sample contributes directly to
#' catalog construction, so systematic per-cohort offsets get absorbed
#' into the compound centroids rather than producing duplicate entries.
#'
#' Each catalog compound gets an extra list-column
#' `per_source_prevalence` recording the fraction of each input
#' cohort's samples in which the compound was observed. This mirrors
#' [merge_catalogs()]'s output for downstream compatibility.
#'
#' @param pp_list Named list of `peak_params` objects.
#' @param eps_rt,eps_cv Catalog clustering eps (fraction of image).
#' @param min_support_frac Support threshold. Applies to the combined
#'   sample set — a compound must span at least this fraction of the
#'   total N samples.
#' @param min_cluster_size Minimum peaks per cluster.
#' @param align Optional cross-cohort drift correction. `FALSE` (default)
#'   skips alignment. `TRUE` calls [align_across_cohorts()] on
#'   `pp_list` with sensible defaults (`eps_rt_loose = 3 * eps_rt`,
#'   `eps_cv_loose = 3 * eps_cv`, `min_support_frac = min_support_frac`,
#'   `iterate = TRUE`, `max_iter = 3L`), then builds the universal
#'   catalog on the aligned cohorts. Corrects BOTH intra-cohort and
#'   inter-cohort drift in one step. Pass a named list to override
#'   any of those. Setting `align = TRUE` is a defensible answer to
#'   reviewer concerns about systematic per-sample or per-batch
#'   drift; the compactness-delta table in `compare_all_catalogs.Rmd`
#'   shows whether it helped on your cohorts.
#' @return A `peak_catalog` object with a `per_source_prevalence`
#'   list-column on `compounds` giving cohort-level prevalence, and
#'   `source_names` / `source_n_samples` fields in `parameters`.
#' @export
build_universal_catalog <- function(pp_list,
                                     eps_rt = 0.01, eps_cv = 0.05,
                                     min_support_frac = 0.10,
                                     min_cluster_size = 2L,
                                     coord_mode = c("fraction", "physical"),
                                     align = FALSE) {
  coord_mode <- match.arg(coord_mode)
  if (!isFALSE(align)) {
    align_args <- if (isTRUE(align)) list() else align
    if (!is.list(align_args))
      stop("build_universal_catalog: `align` must be FALSE, TRUE, or ",
           "a named list of align_across_cohorts parameters.")
    ranges <- function(nm) {
      r <- lapply(pp_list, `[[`, nm)
      if (any(vapply(r, is.null, logical(1)))) NULL else do.call(rbind, r)
    }
    align_args$eps_rt_loose     <- align_args$eps_rt_loose     %||%
      .default_align_eps(eps_rt, ranges("time_ranges"), coord_mode,
                         "build_universal_catalog")
    align_args$eps_cv_loose     <- align_args$eps_cv_loose     %||%
      .default_align_eps(eps_cv, ranges("cv_ranges"), coord_mode,
                         "build_universal_catalog")
    align_args$min_support_frac <- align_args$min_support_frac %||% min_support_frac
    align_args$iterate          <- align_args$iterate          %||% TRUE
    align_args$max_iter         <- align_args$max_iter         %||% 3L
    align_args$verbose          <- align_args$verbose          %||% FALSE
    message("build_universal_catalog: running align_across_cohorts ",
            "(eps_rt_loose = ", align_args$eps_rt_loose,
            ", eps_cv_loose = ", align_args$eps_cv_loose, ") ...")
    pp_list <- do.call(align_across_cohorts,
                        c(list(pp_list = pp_list), align_args))
  }
  # merge_peak_params validates every input and carries all per-peak
  # fields (including rt_seconds_raw / cv_volts_raw) in lockstep.
  merged <- merge_peak_params(pp_list)
  cat <- build_peak_catalog(merged,
    eps_rt = eps_rt, eps_cv = eps_cv,
    min_support_frac = min_support_frac,
    min_cluster_size = min_cluster_size,
    coord_mode = coord_mode)

  # Attach per-cohort prevalence for each compound
  n_cohorts <- length(pp_list)
  cohort_names <- merged$cohort_names
  cohort_n_samples <- vapply(pp_list, function(pp) as.integer(pp$n_samples),
                             integer(1))
  # For each sample_idx in merged, which cohort?
  sample_to_cohort <- merged$sample_cohort

  per_source_prev <- vector("list", nrow(cat$compounds))
  for (i in seq_len(nrow(cat$compounds))) {
    members <- cat$compounds$member_indices[[i]]
    unique_samples <- unique(merged$sample_idx_raw[members])
    cohorts_of_samples <- sample_to_cohort[unique_samples]
    n_per_cohort <- table(factor(cohorts_of_samples,
                                   levels = cohort_names))
    prev_vec <- as.numeric(n_per_cohort) / cohort_n_samples
    names(prev_vec) <- cohort_names
    per_source_prev[[i]] <- prev_vec
  }
  cat$compounds$per_source_prevalence <- per_source_prev

  # Update parameters
  cat$parameters$source_names       <- cohort_names
  cat$parameters$source_n_samples   <- cohort_n_samples
  cat
}

#' Split a catalog into two disjoint sub-catalogs by compound
#'
#' Randomly partitions the catalog's compounds into two sets. Useful
#' for encoder-transfer experiments: pretrain the encoder on
#' synthetic data generated from sub-catalog A, then simulate
#' case/control data from sub-catalog B — biomarker compounds and
#' background compounds in B are all guaranteed to have been unseen
#' during pretraining.
#'
#' Small clusters and singletons are left untouched (attached to the
#' first split only) since they're not typically used for synthetic
#' generation.
#'
#' @param catalog A `peak_catalog` object.
#' @param frac Fraction of compounds to place in the first split.
#'   Default `0.5`.
#' @param seed Optional RNG seed for reproducible partitioning.
#' @return A list with `a` and `b`, each a valid `peak_catalog` with
#'   its own compound subset. `a` also carries any `small_clusters`
#'   and `singletons` from the original.
#' @export
split_catalog <- function(catalog, frac = 0.5, seed = NULL) {
  stopifnot(inherits(catalog, "peak_catalog"))
  if (frac <= 0 || frac >= 1) stop("frac must be in (0, 1).")
  if (!is.null(seed)) set.seed(seed)

  compounds <- catalog$compounds
  n <- nrow(compounds)
  n_a <- max(1L, round(frac * n))
  idx_a <- sort(sample.int(n, size = n_a, replace = FALSE))
  idx_b <- setdiff(seq_len(n), idx_a)

  compounds_a <- compounds[idx_a, ]
  compounds_b <- compounds[idx_b, ]
  compounds_a$compound_id <- seq_len(nrow(compounds_a))
  compounds_b$compound_id <- seq_len(nrow(compounds_b))

  base_shell <- list(
    small_clusters = catalog$small_clusters[0, ],
    singletons     = catalog$singletons[0, ],
    noise          = catalog$noise,
    parameters     = catalog$parameters
  )
  cat_a <- structure(
    c(list(compounds = compounds_a),
      base_shell,
      list(small_clusters = catalog$small_clusters,
           singletons     = catalog$singletons))[
        c("compounds", "small_clusters", "singletons",
          "noise", "parameters")],
    class = "peak_catalog"
  )
  cat_b <- structure(
    list(
      compounds      = compounds_b,
      small_clusters = catalog$small_clusters[0, ],
      singletons     = catalog$singletons[0, ],
      noise          = catalog$noise,
      parameters     = catalog$parameters
    ),
    class = "peak_catalog"
  )
  list(a = cat_a, b = cat_b)
}

#' Plot the spatial distribution of a peak catalog
#'
#' Compound positions are shown as red points sized by prevalence,
#' overlaid on a hexbin density map of the underlying pool peaks.
#' Small clusters and singletons are optionally overlaid as blue
#' points for context.
#'
#' Requires `peak_params` (the object used to build the catalog) so
#' we can draw the underlying peak density.
#'
#' @param catalog A `peak_catalog` object.
#' @param peak_params The `peak_params` used to build the catalog.
#' @param show_small,show_singletons Overlay the below-threshold
#'   clusters and singletons? Default `TRUE` for `small`, `FALSE` for
#'   `singletons` (there are usually many).
#'
#' @return A ggplot object.
#' @export
plot_catalog <- function(catalog, peak_params,
                          show_small = TRUE,
                          show_singletons = FALSE) {
  stopifnot(inherits(catalog, "peak_catalog"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_catalog().")
  }

  pool_df <- tibble::tibble(
    rt_loc = peak_params$rt_loc_raw,
    cv_loc = peak_params$cv_loc_raw
  )
  cmp_df <- catalog$compounds
  sml_df <- catalog$small_clusters
  sng_df <- catalog$singletons

  p <- ggplot2::ggplot() +
    ggplot2::geom_hex(data = pool_df,
                       ggplot2::aes(x = cv_loc, y = rt_loc),
                       bins = c(30, 60), alpha = 0.4) +
    ggplot2::scale_fill_gradient(low = "grey90", high = "grey30",
                                   name = "pool peaks")
  if (show_singletons && nrow(sng_df) > 0L) {
    p <- p + ggplot2::geom_point(data = sng_df,
                                   ggplot2::aes(x = cv_loc, y = rt_loc),
                                   color = "#7570b3", size = 0.6,
                                   alpha = 0.4)
  }
  if (show_small && nrow(sml_df) > 0L) {
    p <- p + ggplot2::geom_point(data = sml_df,
                                   ggplot2::aes(x = cv_loc, y = rt_loc,
                                                 size = n_observations),
                                   color = "#1f77b4", alpha = 0.6)
  }
  if (nrow(cmp_df) > 0L) {
    p <- p + ggplot2::geom_point(data = cmp_df,
                                   ggplot2::aes(x = cv_loc, y = rt_loc,
                                                 size = prevalence),
                                   color = "#B2182B", alpha = 0.85) +
      ggplot2::scale_size_continuous(name = "prevalence",
                                       range = c(0.5, 4))
  }
  p + ggplot2::scale_y_reverse() +
      ggplot2::labs(x = "CV (fraction)", y = "RT (fraction)",
                     title = sprintf(
                       "Peak catalog (compounds: %d, small: %d, singletons: %d)",
                       nrow(cmp_df), nrow(sml_df), nrow(sng_df))) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
}
