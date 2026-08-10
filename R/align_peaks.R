# classydms: sample-level peak alignment
# Correct systematic per-sample shifts by identifying anchor compounds
# and shifting each sample so its peaks align with the anchor positions.
# Improves catalog compactness — otherwise even small instrument drift
# can push peaks into different catalog clusters than they should be.

#' Align samples to anchor peak positions
#'
#' Corrects per-sample systematic shifts in `(rt_loc, cv_loc)`. If no
#' `reference_positions` are provided, a preliminary catalog is built
#' internally and its high-prevalence compounds (>= `min_prevalence_for_anchor`)
#' are used as anchors. For each sample:
#'
#' \enumerate{
#'   \item Match each anchor to the nearest peak in that sample within
#'     `(match_radius_rt, match_radius_cv)`
#'   \item Compute the sample's median (`dr`, `dc`) from those matches
#'   \item Apply `-(dr, dc)` to every peak in that sample
#' }
#'
#' The returned `peak_params` has updated `rt_loc_raw` and `cv_loc_raw`
#' plus an `alignment` field recording per-sample shifts and match
#' counts for diagnostics.
#'
#' Optionally iterates: rebuild the catalog on aligned peaks, extract
#' new anchor positions, re-align. Usually converges in 2 iterations.
#'
#' @param peak_params `peak_params` object with `sample_idx_raw`.
#' @param reference_positions Optional data frame with `rt_loc` and
#'   `cv_loc` columns giving anchor positions. If `NULL`, anchors are
#'   derived from a preliminary catalog build.
#' @param min_prevalence_for_anchor Prevalence threshold on the
#'   preliminary catalog for a compound to serve as an anchor.
#' @param match_radius_rt,match_radius_cv Match distance for
#'   sample-peak to anchor pairing (in fraction of image dims). Should
#'   be roughly the size of expected shifts — too tight misses matches;
#'   too loose picks wrong peaks.
#' @param min_anchors_per_sample Sample won't be shifted if fewer
#'   anchors match. A warning summarizes such samples.
#' @param iterate Whether to re-run alignment after rebuilding the
#'   catalog on aligned peaks. Default `TRUE`.
#' @param max_iter Maximum iterations.
#' @param prelim_catalog_eps_rt,prelim_catalog_eps_cv Catalog eps used
#'   for the internal preliminary build. Must be **larger than the
#'   expected drift magnitude** — otherwise drifted samples' peaks form
#'   their own separate clusters and never make it to anchor status.
#'   Defaults `(0.02, 0.05)` accommodate up to ~2% drift (~28 RT
#'   pixels at native scale), which covers most GC-DMS instrument
#'   variation.
#' @param verbose Whether to print progress.
#' @return Updated `peak_params` with aligned `rt_loc_raw`,
#'   `cv_loc_raw`, and an `alignment` field:
#' \describe{
#'   \item{shifts_rt, shifts_cv}{per-sample shifts applied (fraction).}
#'   \item{n_anchors_matched}{per-sample count of matched anchors.}
#'   \item{n_anchors_total}{total anchors used.}
#'   \item{iterations}{number of iterations run.}
#'   \item{converged}{whether shifts stabilized before max_iter.}
#' }
#' @export
align_peak_params <- function(peak_params,
                                reference_positions = NULL,
                                min_prevalence_for_anchor = 0.75,
                                match_radius_rt = 0.03,
                                match_radius_cv = 0.05,
                                min_anchors_per_sample = 5L,
                                iterate = TRUE,
                                max_iter = 3L,
                                prelim_catalog_eps_rt = 0.02,
                                prelim_catalog_eps_cv = 0.05,
                                verbose = TRUE) {
  needed <- c("rt_loc_raw", "cv_loc_raw", "sample_idx_raw",
              "sigma_rt_raw", "sigma_cv_raw", "intensity_raw", "n_samples")
  missing <- needed[vapply(needed, function(nm) is.null(peak_params[[nm]]),
                            logical(1))]
  if (length(missing)) {
    stop("align_peak_params: peak_params is missing ",
         paste(missing, collapse = ", "),
         ". Call backfill_sample_idx_raw() if needed.")
  }
  n_samples <- peak_params$n_samples

  # Track total accumulated shifts across iterations
  total_shift_rt <- numeric(n_samples)
  total_shift_cv <- numeric(n_samples)
  converged <- FALSE
  iterations_run <- 0L

  current_pp <- peak_params
  prev_iter_shifts <- NULL

  for (iter in seq_len(max_iter)) {
    iterations_run <- iter

    # Get anchors: from reference_positions on first iteration if given,
    # else from preliminary catalog on current (possibly-aligned) peaks
    anchors <- if (iter == 1L && !is.null(reference_positions)) {
      reference_positions
    } else {
      pc <- build_peak_catalog(current_pp,
        eps_rt = prelim_catalog_eps_rt,
        eps_cv = prelim_catalog_eps_cv,
        min_support_frac = 0.10)
      anchors_df <- pc$compounds[
        pc$compounds$prevalence >= min_prevalence_for_anchor,
        c("rt_loc", "cv_loc")]
      if (nrow(anchors_df) < min_anchors_per_sample) {
        warning("Only ", nrow(anchors_df),
                " anchor candidates at prevalence >= ",
                min_prevalence_for_anchor,
                ". Try lowering min_prevalence_for_anchor.")
      }
      anchors_df
    }
    if (verbose) {
      message(sprintf("  Iter %d: using %d anchors", iter, nrow(anchors)))
    }

    rt <- current_pp$rt_loc_raw
    cv <- current_pp$cv_loc_raw
    sample_idx <- current_pp$sample_idx_raw
    aligned_rt <- rt
    aligned_cv <- cv
    shifts_rt <- numeric(n_samples)
    shifts_cv <- numeric(n_samples)
    n_matched <- integer(n_samples)

    for (s in seq_len(n_samples)) {
      peak_mask <- sample_idx == s
      if (!any(peak_mask)) next
      peaks_rt <- rt[peak_mask]
      peaks_cv <- cv[peak_mask]

      # For each anchor, collect ALL sample peaks within match radius —
      # not just the nearest. The true drift-shifted peak forms a
      # consistent pile at (dr = drift), while random unrelated peaks
      # scatter uniformly across the window. Median across all pairs
      # collapses onto the mode where the pile is.
      dr_list <- numeric(0)
      dc_list <- numeric(0)
      n_anchors_hit <- 0L
      for (a in seq_len(nrow(anchors))) {
        dr <- peaks_rt - anchors$rt_loc[a]
        dc <- peaks_cv - anchors$cv_loc[a]
        within <- which(abs(dr) <= match_radius_rt &
                         abs(dc) <= match_radius_cv)
        if (length(within) > 0L) {
          dr_list <- c(dr_list, dr[within])
          dc_list <- c(dc_list, dc[within])
          n_anchors_hit <- n_anchors_hit + 1L
        }
      }

      n_matched[s] <- n_anchors_hit
      if (n_anchors_hit >= min_anchors_per_sample &&
          length(dr_list) > 0L) {
        shifts_rt[s] <- stats::median(dr_list)
        shifts_cv[s] <- stats::median(dc_list)
        aligned_rt[peak_mask] <- peaks_rt - shifts_rt[s]
        aligned_cv[peak_mask] <- peaks_cv - shifts_cv[s]
      }
    }

    # Update running total shifts (which reflect the ORIGINAL positions)
    total_shift_rt <- total_shift_rt + shifts_rt
    total_shift_cv <- total_shift_cv + shifts_cv

    current_pp$rt_loc_raw <- aligned_rt
    current_pp$cv_loc_raw <- aligned_cv

    # Check convergence: max |shift| this iteration
    max_delta <- max(abs(shifts_rt), abs(shifts_cv))
    if (verbose) {
      message(sprintf(
        "    max |shift| this iter: %.5f  (median RT: %.5f, CV: %.5f)",
        max_delta, stats::median(shifts_rt), stats::median(shifts_cv)))
    }
    if (!iterate) break
    # Consider converged if max shift is below 1/10 of match radius
    if (max_delta < min(match_radius_rt, match_radius_cv) / 10) {
      converged <- TRUE
      if (verbose) message("    converged")
      break
    }
    prev_iter_shifts <- shifts_rt
  }

  # Update derived stats
  current_pp$rt_loc <- list(
    values = current_pp$rt_loc_raw,
    mean = mean(current_pp$rt_loc_raw),
    sd = stats::sd(current_pp$rt_loc_raw))
  current_pp$cv_loc <- list(
    values = current_pp$cv_loc_raw,
    mean = mean(current_pp$cv_loc_raw),
    sd = stats::sd(current_pp$cv_loc_raw))

  current_pp$alignment <- list(
    shifts_rt = total_shift_rt,
    shifts_cv = total_shift_cv,
    n_anchors_matched = n_matched,
    n_anchors_total   = nrow(anchors),
    iterations        = iterations_run,
    converged         = converged
  )

  if (verbose) {
    cat(sprintf(
      "align_peak_params: %d samples aligned via %d anchors, %d iters\n",
      n_samples, nrow(anchors), iterations_run))
    cat(sprintf(
      "  total shift: median RT = %.4f, CV = %.4f\n",
      stats::median(total_shift_rt), stats::median(total_shift_cv)))
    cat(sprintf(
      "  max |total shift|: RT = %.4f, CV = %.4f\n",
      max(abs(total_shift_rt)), max(abs(total_shift_cv))))
    under <- sum(n_matched < min_anchors_per_sample)
    if (under > 0) {
      cat(sprintf("  WARNING: %d samples had < %d anchor matches (not shifted)\n",
                  under, min_anchors_per_sample))
    }
  }
  current_pp
}
