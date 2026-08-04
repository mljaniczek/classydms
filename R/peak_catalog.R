# classydms: peak catalog construction and diagnostics
# Groups pool peaks into "compounds" — specific (RT, CV) positions
# where peaks reliably appear across samples — for realistic synthetic
# data generation. See PIPELINE_KNOBS.md and the case-control tutorial
# for context.

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
#'     intensity_obs, member_indices`.}
#'   \item{small_clusters}{same schema, for clusters below thresholds.}
#'   \item{singletons}{same schema, for clusters of size 1.}
#'   \item{parameters}{the settings used to build the catalog.}
#' }
#' @export
build_peak_catalog <- function(peak_params,
                                eps_rt = 0.01, eps_cv = 0.05,
                                min_support_frac = 0.10,
                                min_cluster_size = 2L) {
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

  rt         <- peak_params$rt_loc_raw
  cv         <- peak_params$cv_loc_raw
  intensity  <- peak_params$intensity_raw
  sample_idx <- peak_params$sample_idx_raw
  n_peaks    <- length(rt)
  n_samples  <- peak_params$n_samples

  if (n_peaks == 0L) {
    stop("build_peak_catalog: peak_params$rt_loc_raw is empty.")
  }

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

  cluster_rows <- lapply(seq_len(n_clusters), function(cid) {
    members <- cluster_members[[as.character(cid)]]
    unique_samples <- unique(sample_idx[members])
    tibble::tibble(
      compound_id     = cid,
      rt_loc          = stats::median(rt[members]),
      cv_loc          = stats::median(cv[members]),
      location_sd_rt  = if (length(members) > 1L)
                          stats::sd(rt[members]) else 0,
      location_sd_cv  = if (length(members) > 1L)
                          stats::sd(cv[members]) else 0,
      n_observations  = length(members),
      n_samples       = length(unique_samples),
      prevalence      = length(unique_samples) / n_samples,
      sigma_rt_obs    = list(peak_params$sigma_rt_raw[members]),
      sigma_cv_obs    = list(peak_params$sigma_cv_raw[members]),
      intensity_obs   = list(intensity[members]),
      member_indices  = list(members)
    )
  })
  clusters_df <- dplyr::bind_rows(cluster_rows)

  # Classify into compounds / small_clusters / singletons.
  is_singleton      <- clusters_df$n_observations < 2L
  is_below_support  <- clusters_df$prevalence < min_support_frac
  is_below_size     <- clusters_df$n_observations < min_cluster_size
  is_compound       <- !is_singleton & !is_below_support & !is_below_size

  compounds       <- clusters_df[is_compound, ]
  small_clusters  <- clusters_df[!is_compound & !is_singleton, ]
  singletons      <- clusters_df[is_singleton, ]

  # Re-number compounds starting from 1 for convenience.
  if (nrow(compounds) > 0L) {
    compounds$compound_id <- seq_len(nrow(compounds))
  }

  structure(
    list(
      compounds       = compounds,
      small_clusters  = small_clusters,
      singletons      = singletons,
      parameters      = list(
        eps_rt              = eps_rt,
        eps_cv              = eps_cv,
        min_support_frac    = min_support_frac,
        min_cluster_size    = min_cluster_size,
        n_pool_peaks        = n_peaks,
        n_samples           = n_samples,
        n_total_clusters    = n_clusters
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
