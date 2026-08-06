# classydms: cross-catalog operations
# Compare and merge peak catalogs derived from different subsets
# (e.g. per-disease, per-batch, disease vs control).

#' Compare two peak catalogs by compound proximity
#'
#' Matches each compound in `catalog_a` to its nearest neighbor in
#' `catalog_b` by spatial proximity in fraction-of-image `(rt_loc,
#' cv_loc)` space. Compounds are "matched" if the nearest neighbor is
#' within `(match_radius_rt, match_radius_cv)`. Returns three data
#' frames: shared compounds (matched to each other), compounds only in
#' A, and compounds only in B. Also reports the per-catalog prevalence
#' for shared compounds so you can see whether the same compound
#' expresses at different rates between cohorts.
#'
#' @param catalog_a,catalog_b `peak_catalog` objects.
#' @param name_a,name_b Short labels used as column-name suffixes in
#'   the `shared` output. Default `"a"` and `"b"`.
#' @param match_radius_rt,match_radius_cv Maximum distance to count a
#'   compound as matched. Defaults `(0.01, 0.05)` — the same defaults
#'   as `build_peak_catalog`.
#'
#' @return List with:
#' \describe{
#'   \item{shared}{tibble of matched compound pairs with
#'     `compound_id_<a>`, `compound_id_<b>`, `rt_loc_<a>`,
#'     `rt_loc_<b>`, `cv_loc_<a>`, `cv_loc_<b>`,
#'     `prevalence_<a>`, `prevalence_<b>`, `prevalence_delta`
#'     (b − a), and `spatial_distance` (Euclidean in the fraction
#'     space).}
#'   \item{only_a}{tibble of compounds present in A but not matched
#'     in B.}
#'   \item{only_b}{tibble of compounds present in B but not matched
#'     in A.}
#'   \item{summary}{one-row tibble of counts + fractions.}
#' }
#' @export
compare_catalogs <- function(catalog_a, catalog_b,
                              name_a = "a", name_b = "b",
                              match_radius_rt = 0.01,
                              match_radius_cv = 0.05) {
  stopifnot(inherits(catalog_a, "peak_catalog"),
             inherits(catalog_b, "peak_catalog"))
  ca <- catalog_a$compounds
  cb <- catalog_b$compounds
  if (nrow(ca) == 0L || nrow(cb) == 0L) {
    stop("compare_catalogs: at least one catalog has zero compounds.")
  }

  # For each compound in a, find the nearest compound in b within radii.
  matched_b_idx  <- integer(nrow(ca))
  matched_b_dist <- rep(NA_real_, nrow(ca))
  for (i in seq_len(nrow(ca))) {
    dr <- abs(cb$rt_loc - ca$rt_loc[i])
    dc <- abs(cb$cv_loc - ca$cv_loc[i])
    within <- which(dr <= match_radius_rt & dc <= match_radius_cv)
    if (length(within) == 0L) {
      matched_b_idx[i]  <- NA_integer_
      next
    }
    # Nearest by Euclidean in normalized fraction space.
    d <- sqrt((dr[within] / match_radius_rt)^2 +
              (dc[within] / match_radius_cv)^2)
    j <- within[which.min(d)]
    matched_b_idx[i]  <- j
    matched_b_dist[i] <- sqrt(dr[j]^2 + dc[j]^2)
  }

  matched_a <- which(!is.na(matched_b_idx))
  matched_b <- matched_b_idx[matched_a]

  make_col_name <- function(base, suffix) paste0(base, "_", suffix)
  shared <- tibble::tibble(
    !!make_col_name("compound_id", name_a) := ca$compound_id[matched_a],
    !!make_col_name("compound_id", name_b) := cb$compound_id[matched_b],
    !!make_col_name("rt_loc", name_a)      := ca$rt_loc[matched_a],
    !!make_col_name("rt_loc", name_b)      := cb$rt_loc[matched_b],
    !!make_col_name("cv_loc", name_a)      := ca$cv_loc[matched_a],
    !!make_col_name("cv_loc", name_b)      := cb$cv_loc[matched_b],
    !!make_col_name("prevalence", name_a)  := ca$prevalence[matched_a],
    !!make_col_name("prevalence", name_b)  := cb$prevalence[matched_b],
    prevalence_delta = cb$prevalence[matched_b] - ca$prevalence[matched_a],
    spatial_distance = matched_b_dist[matched_a]
  )

  only_a <- ca[is.na(matched_b_idx),
                c("compound_id", "rt_loc", "cv_loc",
                  "prevalence", "n_observations", "n_samples")]
  b_matched_set <- unique(matched_b)
  only_b <- cb[!(seq_len(nrow(cb)) %in% b_matched_set),
                c("compound_id", "rt_loc", "cv_loc",
                  "prevalence", "n_observations", "n_samples")]

  summary <- tibble::tibble(
    n_a          = nrow(ca),
    n_b          = nrow(cb),
    n_shared     = nrow(shared),
    n_only_a     = nrow(only_a),
    n_only_b     = nrow(only_b),
    frac_shared_a = nrow(shared) / nrow(ca),
    frac_shared_b = nrow(shared) / nrow(cb)
  )

  list(shared = shared, only_a = only_a, only_b = only_b,
        summary = summary)
}

#' Merge multiple peak catalogs into a single unified catalog
#'
#' Combines compounds across catalogs by spatial proximity. Overlapping
#' compounds (within `(match_radius_rt, match_radius_cv)` of each
#' other across catalogs) are merged into a single entry: observations
#' are concatenated, and prevalence is recomputed as a weighted average
#' by cohort size. Compounds unique to a single catalog are added as
#' new entries.
#'
#' The merged catalog has an extra `per_source_prevalence` list-column
#' recording per-input-catalog prevalence for each merged compound —
#' useful for tracking which sources contributed which compounds.
#'
#' @param catalog_list List of `peak_catalog` objects.
#' @param names Optional character vector of length `length(catalog_list)`
#'   naming each source. Default `"catalog_1"`, `"catalog_2"`, etc.
#' @param match_radius_rt,match_radius_cv Cross-catalog matching radii.
#'   Defaults to the minimum `eps_rt` / `eps_cv` across the source
#'   catalogs — using a larger radius than the source eps would cause
#'   within-source compounds to fold into each other during the merge.
#'
#' @return A `peak_catalog` object with the merged compounds.
#' @export
merge_catalogs <- function(catalog_list,
                            names = NULL,
                            match_radius_rt = NULL,
                            match_radius_cv = NULL) {
  stopifnot(is.list(catalog_list))
  if (!all(vapply(catalog_list, inherits, logical(1), "peak_catalog"))) {
    stop("merge_catalogs: all elements must be peak_catalog objects.")
  }
  n_cats <- length(catalog_list)
  if (is.null(names))
    names <- paste0("catalog_", seq_len(n_cats))
  stopifnot(length(names) == n_cats)

  # Default match radii to min(source eps) so within-source compounds
  # can't fold into each other during merge.
  eps_rts <- vapply(catalog_list,
                     function(c) c$parameters$eps_rt, numeric(1))
  eps_cvs <- vapply(catalog_list,
                     function(c) c$parameters$eps_cv, numeric(1))
  if (is.null(match_radius_rt)) match_radius_rt <- min(eps_rts)
  if (is.null(match_radius_cv)) match_radius_cv <- min(eps_cvs)

  # Grab per-catalog sample counts for prevalence weighting
  n_samples_vec <- vapply(catalog_list,
                           function(c) c$parameters$n_samples,
                           integer(1))

  # Accumulator: growing list of merged compounds with per-source
  # prevalence contributions.
  merged_rt      <- numeric(0)
  merged_cv      <- numeric(0)
  merged_sigma_rt_obs   <- list()
  merged_sigma_cv_obs   <- list()
  merged_intensity_obs  <- list()
  merged_per_source_prev <- list()   # each element is a length-n_cats vector
  merged_n_obs   <- integer(0)

  for (k in seq_len(n_cats)) {
    cat_k <- catalog_list[[k]]$compounds
    for (i in seq_len(nrow(cat_k))) {
      rt_i <- cat_k$rt_loc[i]; cv_i <- cat_k$cv_loc[i]
      # Look for a match in the growing merged set
      if (length(merged_rt) > 0L) {
        dr <- abs(merged_rt - rt_i); dc <- abs(merged_cv - cv_i)
        within <- which(dr <= match_radius_rt & dc <= match_radius_cv)
      } else {
        within <- integer(0)
      }
      if (length(within) > 0L) {
        # Attach to closest merged compound.
        d <- sqrt((dr[within] / match_radius_rt)^2 +
                  (dc[within] / match_radius_cv)^2)
        m <- within[which.min(d)]
        merged_sigma_rt_obs[[m]]  <- c(merged_sigma_rt_obs[[m]],
                                        cat_k$sigma_rt_obs[[i]])
        merged_sigma_cv_obs[[m]]  <- c(merged_sigma_cv_obs[[m]],
                                        cat_k$sigma_cv_obs[[i]])
        merged_intensity_obs[[m]] <- c(merged_intensity_obs[[m]],
                                        cat_k$intensity_obs[[i]])
        merged_per_source_prev[[m]][k] <- cat_k$prevalence[i]
        merged_n_obs[m] <- merged_n_obs[m] + cat_k$n_observations[i]
        # Update centroid to weighted mean
        w1 <- merged_n_obs[m] - cat_k$n_observations[i]
        w2 <- cat_k$n_observations[i]
        merged_rt[m] <- (merged_rt[m] * w1 + rt_i * w2) / (w1 + w2)
        merged_cv[m] <- (merged_cv[m] * w1 + cv_i * w2) / (w1 + w2)
      } else {
        # New compound
        merged_rt      <- c(merged_rt, rt_i)
        merged_cv      <- c(merged_cv, cv_i)
        merged_sigma_rt_obs[[length(merged_rt)]]  <- cat_k$sigma_rt_obs[[i]]
        merged_sigma_cv_obs[[length(merged_rt)]]  <- cat_k$sigma_cv_obs[[i]]
        merged_intensity_obs[[length(merged_rt)]] <- cat_k$intensity_obs[[i]]
        prev_vec <- rep(0, n_cats)
        prev_vec[k] <- cat_k$prevalence[i]
        merged_per_source_prev[[length(merged_rt)]] <- prev_vec
        merged_n_obs   <- c(merged_n_obs, cat_k$n_observations[i])
      }
    }
  }

  # Combined prevalence = weighted mean of per-source prevalences
  # (weighted by n_samples per source).
  total_samples <- sum(n_samples_vec)
  combined_prev <- vapply(merged_per_source_prev, function(v)
    sum(v * n_samples_vec) / total_samples,
    numeric(1))
  n_samples_total <- vapply(merged_per_source_prev, function(v)
    sum((v > 0) * n_samples_vec),
    numeric(1))

  merged_df <- tibble::tibble(
    compound_id     = seq_along(merged_rt),
    rt_loc          = merged_rt,
    cv_loc          = merged_cv,
    location_sd_rt  = 0,   # not meaningful across catalogs
    location_sd_cv  = 0,
    n_observations  = merged_n_obs,
    n_samples       = as.integer(n_samples_total),
    prevalence      = combined_prev,
    sigma_rt_obs    = merged_sigma_rt_obs,
    sigma_cv_obs    = merged_sigma_cv_obs,
    intensity_obs   = merged_intensity_obs,
    member_indices  = replicate(length(merged_rt),
                                  integer(0), simplify = FALSE),
    per_source_prevalence = lapply(merged_per_source_prev,
      function(v) stats::setNames(v, names))
  )

  # Combine noise stats: mean of per-catalog noise means, mean of SDs.
  noise_means <- vapply(catalog_list,
                         function(c) c$noise$mean, numeric(1))
  noise_sds   <- vapply(catalog_list,
                         function(c) c$noise$sd, numeric(1))
  noise <- list(mean = mean(noise_means, na.rm = TRUE),
                sd   = mean(noise_sds,   na.rm = TRUE))

  structure(
    list(
      compounds       = merged_df,
      small_clusters  = catalog_list[[1]]$small_clusters[0, ],
      singletons      = catalog_list[[1]]$singletons[0, ],
      noise           = noise,
      parameters      = list(
        eps_rt              = match_radius_rt,
        eps_cv              = match_radius_cv,
        min_support_frac    = 0,
        min_cluster_size    = 1L,
        n_pool_peaks        = sum(vapply(catalog_list,
                                     function(c) c$parameters$n_pool_peaks,
                                     integer(1))),
        n_samples           = total_samples,
        n_total_clusters    = length(merged_rt),
        source_names        = names,
        source_n_samples    = n_samples_vec
      )
    ),
    class = "peak_catalog"
  )
}
