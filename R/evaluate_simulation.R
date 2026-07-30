# classydms: evaluation metrics for the case / control simulation
# Given a saliency map produced by the pipeline and a ground-truth
# biomarker spec, report how well saliency localizes the injected
# signal.

#' Non-max-suppressed local maxima of a matrix, ranked by value
#'
#' Very light wrapper on the same greedy NMS used for peak detection —
#' pick the brightest pixel, then the next-brightest that is more than
#' `min_sep_*` away from every accepted pixel on at least one axis, up
#' to `top_k` picks.
#'
#' @param M numeric matrix (typically a saliency map).
#' @param top_k number of local maxima to return.
#' @param min_sep_rt,min_sep_cv separation on each axis (in pixels).
#' @return integer matrix with columns `r`, `c` and `val`, ranked by `val`.
#' @keywords internal
top_local_maxima <- function(M, top_k = 10L,
                              min_sep_rt = 8L, min_sep_cv = 2L) {
  H <- nrow(M); W <- ncol(M)
  idx <- which(M > 0)
  if (length(idx) == 0L) {
    return(data.frame(r = integer(0), c = integer(0), val = numeric(0)))
  }
  ord <- idx[order(M[idx], decreasing = TRUE)]
  centers <- data.frame(r = integer(0), c = integer(0), val = numeric(0))
  for (k in seq_along(ord)) {
    if (nrow(centers) >= top_k) break
    rc <- arrayInd(ord[k], .dim = c(H, W))
    r0 <- rc[1]; c0 <- rc[2]
    if (nrow(centers) == 0L) {
      centers <- rbind(centers, data.frame(r = r0, c = c0, val = M[r0, c0]))
      next
    }
    dr <- abs(centers$r - r0); dc <- abs(centers$c - c0)
    if (all(dr >= min_sep_rt | dc >= min_sep_cv)) {
      centers <- rbind(centers, data.frame(r = r0, c = c0, val = M[r0, c0]))
    }
  }
  centers
}

#' Evaluate biomarker recovery from a saliency map
#'
#' Two complementary metrics on the same saliency map and ground-truth
#' biomarker spec:
#'
#' \enumerate{
#'   \item \strong{Per-pixel AUC}: label each pixel `1` if within
#'     `(match_radius_rt, match_radius_cv)` of any biomarker center,
#'     else `0`; score with the saliency value; report ROC AUC. This
#'     is "does saliency track biomarker regions on average?"
#'   \item \strong{Top-K precision / recall}: extract the top-K
#'     non-max-suppressed local maxima of the saliency map; for each,
#'     check whether it falls within the match radius of any biomarker
#'     center. Precision = matched / K. Recall = unique biomarkers
#'     matched / n_biomarkers. This is "does saliency's hotspot list
#'     correspond to the biomarkers?"
#' }
#'
#' @param saliency_map numeric `H` x `W` matrix in the same coordinate
#'   frame as the biomarker locations (i.e. after any padding /
#'   trimming used for the encoder input).
#' @param biomarkers ground-truth spec: tibble with at least `rt_loc`
#'   and `cv_loc` columns (as fractions of image dimensions).
#' @param H,W dimensions of `saliency_map`. If `NULL`, taken from the
#'   matrix.
#' @param match_radius_rt,match_radius_cv match radius in pixels.
#'   Default 8 RT / 2 CV, matching the peak-detector's default
#'   `min_sep_rt`, `min_sep_cv`.
#' @param top_k number of local maxima to extract for precision/recall.
#'
#' @return list(auc_per_pixel, precision_at_k, recall_at_k,
#'   matched_biomarkers, missed_biomarkers, top_hotspots)
#' @export
evaluate_biomarker_recovery <- function(saliency_map,
                                         biomarkers,
                                         H = NULL, W = NULL,
                                         match_radius_rt = 8L,
                                         match_radius_cv = 2L,
                                         top_k = 10L) {
  if (is.null(H)) H <- nrow(saliency_map)
  if (is.null(W)) W <- ncol(saliency_map)
  stopifnot(nrow(saliency_map) == H, ncol(saliency_map) == W)

  bm_r <- biomarkers$rt_loc * H
  bm_c <- biomarkers$cv_loc * W
  n_bm <- length(bm_r)

  # -------- (1) per-pixel AUC ---------------------------------------
  # Ground-truth mask: 1 within rectangle around any biomarker.
  gt_mask <- matrix(0L, nrow = H, ncol = W)
  for (b in seq_len(n_bm)) {
    r_lo <- max(1L, floor(bm_r[b] - match_radius_rt))
    r_hi <- min(H,  ceiling(bm_r[b] + match_radius_rt))
    c_lo <- max(1L, floor(bm_c[b] - match_radius_cv))
    c_hi <- min(W,  ceiling(bm_c[b] + match_radius_cv))
    gt_mask[r_lo:r_hi, c_lo:c_hi] <- 1L
  }
  labels <- as.vector(gt_mask)
  scores <- as.vector(saliency_map)
  auc_per_pixel <- if (length(unique(labels)) < 2L) NA_real_ else {
    tryCatch(
      as.numeric(pROC::auc(pROC::roc(labels, scores,
                                       quiet = TRUE, direction = "<"))),
      error = function(e) NA_real_
    )
  }

  # -------- (2) top-K precision / recall ----------------------------
  hotspots <- top_local_maxima(saliency_map, top_k = top_k,
                                min_sep_rt = match_radius_rt,
                                min_sep_cv = match_radius_cv)

  matched_biomarkers <- integer(0)
  hotspot_matches <- logical(nrow(hotspots))
  if (nrow(hotspots) > 0L) {
    for (h in seq_len(nrow(hotspots))) {
      dr <- abs(bm_r - hotspots$r[h])
      dc <- abs(bm_c - hotspots$c[h])
      hits <- which(dr <= match_radius_rt & dc <= match_radius_cv)
      if (length(hits) > 0L) {
        hotspot_matches[h] <- TRUE
        matched_biomarkers <- union(matched_biomarkers, hits)
      }
    }
  }
  precision_at_k <- if (nrow(hotspots) == 0L) NA_real_
                    else mean(hotspot_matches)
  recall_at_k    <- length(matched_biomarkers) / n_bm
  missed_biomarkers <- setdiff(seq_len(n_bm), matched_biomarkers)

  list(
    auc_per_pixel      = auc_per_pixel,
    precision_at_k     = precision_at_k,
    recall_at_k        = recall_at_k,
    matched_biomarkers = matched_biomarkers,
    missed_biomarkers  = missed_biomarkers,
    top_hotspots       = hotspots,
    gt_mask            = gt_mask
  )
}
