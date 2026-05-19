# classydms: per-sample quality control
# Compute per-sample QC features and flag samples that should be
# excluded from downstream analysis. Designed to be called between
# preprocessing (process_one_sample + baseline_basement + trim_sample)
# and peak parameter estimation / encoder pretraining.

#' Compute per-sample QC features
#'
#' Calculates a battery of sample-level summary statistics on a list of
#' preprocessed sample matrices (typically after Savitzky-Golay + AsLS +
#' dust thresholding + trim, but before normalization or padding).
#' Returns a tibble with one row per sample and columns covering signal
#' intensity, sparsity, peak density, and (optional) cohort-relative
#' similarity.
#'
#' The returned table is consumed by [flag_junk_samples()], which
#' applies user-tunable thresholds to mark samples for exclusion.
#'
#' @param samples A list of sample objects, each with element `Z` (the
#'   preprocessed intensity matrix). Typically the output of
#'   `lapply(samples, trim_sample, bounds = global_b)`.
#' @param sample_names Optional character vector of sample names. If
#'   `NULL`, uses `names(samples)` or falls back to integer indices.
#' @param compute_cohort_correlation If `TRUE`, also compute each
#'   sample's correlation against the cohort median spectrum. Slower
#'   (O(n_samples * n_pixels)) but useful for detecting morphological
#'   outliers that don't trigger univariate flags. Default `FALSE`.
#' @return A tibble with one row per sample and columns:
#'   - `sample` — sample identifier
#'   - `rt_dim`, `cv_dim` — post-trim dimensions
#'   - `n_positive`, `prop_nonzero` — count and fraction of pixels above zero
#'   - `total_signal` — sum of all intensities
#'   - `mean_positive` — mean of pixels above zero
#'   - `max_intensity` — peak pixel value
#'   - `dynamic_range` — `max_intensity / median(positive pixels)`
#'   - `bbox_size` — `rt_dim * cv_dim`
#'   - `n_peaks` — detected peaks via [pick_peak_centers()]
#'   - `cohort_correlation` — (if requested) correlation with cohort median
#' @export
qc_features <- function(samples,
                         sample_names = NULL,
                         compute_cohort_correlation = FALSE) {

  if (is.null(sample_names)) {
    sample_names <- names(samples) %||% as.character(seq_along(samples))
  }

  per_sample <- lapply(seq_along(samples), function(i) {
    Z <- samples[[i]]$Z
    if (is.null(Z)) Z <- samples[[i]]   # support raw matrix list too
    pos <- Z[Z > 0]
    n_pos <- length(pos)

    list(
      sample        = sample_names[i],
      rt_dim        = nrow(Z),
      cv_dim        = ncol(Z),
      n_positive    = n_pos,
      prop_nonzero  = n_pos / length(Z),
      total_signal  = sum(Z),
      mean_positive = if (n_pos > 0) mean(pos) else NA_real_,
      max_intensity = max(Z),
      dynamic_range = if (n_pos > 0) max(Z) / stats::median(pos) else NA_real_,
      bbox_size     = nrow(Z) * ncol(Z),
      n_peaks = {
        if (n_pos > 0) {
          eps <- 0.01 * stats::quantile(pos, 0.95)
          nrow(pick_peak_centers(Z, top_k = 10000L, eps = eps,
                                  min_sep_rt = 8L, min_sep_cv = 2L))
        } else {
          0L
        }
      }
    )
  })

  qc <- do.call(rbind, lapply(per_sample, as.data.frame,
                               stringsAsFactors = FALSE))
  qc <- tibble::as_tibble(qc)

  if (compute_cohort_correlation) {
    # Resample each sample to a common low-res grid, then compute
    # correlation of each sample's flattened spectrum against the cohort
    # median. The downsampling is purely to make this O(n_samples * 1000)
    # rather than O(n_samples * 100000).
    target_dim <- c(64L, 32L)
    flatten_one <- function(Z) {
      # Simple block average to target dim
      r_idx <- pmin(ceiling(seq(1, nrow(Z), length.out = target_dim[1])),
                    nrow(Z))
      c_idx <- pmin(ceiling(seq(1, ncol(Z), length.out = target_dim[2])),
                    ncol(Z))
      as.numeric(Z[r_idx, c_idx])
    }
    flat_mat <- sapply(samples, function(s) {
      Z <- s$Z; if (is.null(Z)) Z <- s
      flatten_one(Z)
    })  # [pixels x samples]
    cohort_median <- apply(flat_mat, 1, stats::median, na.rm = TRUE)
    qc$cohort_correlation <- vapply(seq_len(ncol(flat_mat)), function(i) {
      v <- flat_mat[, i]
      if (stats::sd(v) == 0 || stats::sd(cohort_median) == 0) return(NA_real_)
      stats::cor(v, cohort_median, use = "complete.obs")
    }, numeric(1))
  }

  qc
}

#' Flag samples that fail QC checks
#'
#' Applies threshold-based and IQR-based criteria to a QC table from
#' [qc_features()] and returns the table with boolean flag columns
#' added. A sample is marked `flagged = TRUE` if it triggers any of
#' the active criteria. Thresholds are tunable; pass `NA` or `Inf` to
#' disable a specific check.
#'
#' Default criteria (Tier 1 + Tier 2 from the package design):
#'
#' - **Low signal**: `prop_nonzero` below `min_prop_nonzero`.
#' - **Saturation**: `max_intensity` at or above `max_saturation`.
#' - **Low dynamic range**: `dynamic_range` below `min_dynamic_range`.
#' - **Univariate IQR outliers** on `total_signal`, `mean_positive`,
#'   `n_peaks`, `prop_nonzero`, `max_intensity`, `dynamic_range`, and
#'   `bbox_size` (any feature more than `iqr_mult` IQR from the cohort
#'   median triggers).
#' - **Low cohort correlation** (only if `cohort_correlation` is in
#'   the table): correlation below `min_cohort_correlation`.
#'
#' @param qc_tbl Output of [qc_features()].
#' @param min_prop_nonzero Threshold below which a sample is flagged as
#'   low-signal (default 0.10).
#' @param max_saturation Threshold at or above which a sample is flagged
#'   for saturation (default `Inf` = disabled; set to e.g. 65535 for
#'   16-bit detectors, or whatever the instrument's max readout is).
#' @param min_dynamic_range Threshold below which dynamic range is
#'   considered too low (default 5).
#' @param iqr_mult Multiplier for IQR-based outlier detection across
#'   univariate features (default 3).
#' @param min_cohort_correlation Threshold below which a sample's
#'   correlation with the cohort median spectrum triggers a flag
#'   (default 0.3; only applied if `cohort_correlation` is in the table).
#' @param features_for_iqr Character vector of QC feature names to
#'   include in IQR outlier detection. Default uses the package-recommended
#'   set; pass NULL to disable IQR checks entirely.
#' @return The input tibble with additional columns: one `*_outlier`
#'   column per IQR-checked feature, plus `low_signal`, `saturation`,
#'   `low_dynamic_range`, `low_cohort_correlation`, and a final
#'   `flagged` column that is the OR of all criteria.
#' @export
flag_junk_samples <- function(qc_tbl,
                                min_prop_nonzero = 0.10,
                                max_saturation = Inf,
                                min_dynamic_range = 5,
                                iqr_mult = 3,
                                min_cohort_correlation = 0.30,
                                features_for_iqr = c("total_signal",
                                                      "mean_positive",
                                                      "n_peaks",
                                                      "prop_nonzero",
                                                      "max_intensity",
                                                      "dynamic_range",
                                                      "bbox_size")) {

  flags <- qc_tbl
  flags$low_signal       <- flags$prop_nonzero < min_prop_nonzero
  flags$saturation       <- is.finite(max_saturation) &
                              flags$max_intensity >= max_saturation
  flags$low_dynamic_range <- !is.na(flags$dynamic_range) &
                              flags$dynamic_range < min_dynamic_range

  if (!is.null(features_for_iqr)) {
    for (col in features_for_iqr) {
      if (!col %in% names(flags)) next
      v <- flags[[col]]
      q1 <- stats::quantile(v, 0.25, na.rm = TRUE)
      q3 <- stats::quantile(v, 0.75, na.rm = TRUE)
      iqr <- q3 - q1
      flags[[paste0(col, "_outlier")]] <-
        !is.na(v) &
        (v < q1 - iqr_mult * iqr | v > q3 + iqr_mult * iqr)
    }
  }

  if ("cohort_correlation" %in% names(flags)) {
    flags$low_cohort_correlation <- !is.na(flags$cohort_correlation) &
                                      flags$cohort_correlation <
                                        min_cohort_correlation
  }

  # Combine all flag columns into the final `flagged` indicator
  flag_cols <- grep(
    "^(low_signal|saturation|low_dynamic_range|low_cohort_correlation|.*_outlier)$",
    names(flags), value = TRUE)
  flags$flagged <- Reduce("|", lapply(flag_cols, function(c) flags[[c]]))

  flags
}

#' Visualize per-sample QC features
#'
#' Convenience plot showing the cohort distribution of each QC feature,
#' with flagged samples highlighted. Useful for inspecting whether the
#' default thresholds are reasonable for your cohort before committing
#' to exclusion.
#'
#' @param qc_flagged Output of [flag_junk_samples()].
#' @param features Character vector of features to plot (default: the
#'   IQR-checked set).
#' @return A ggplot object (one panel per feature).
#' @export
plot_qc <- function(qc_flagged,
                     features = c("total_signal", "mean_positive",
                                  "n_peaks", "prop_nonzero",
                                  "max_intensity", "dynamic_range")) {

  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE)) {
    stop("Requires ggplot2 and tidyr (in Suggests).")
  }

  available <- features[features %in% names(qc_flagged)]
  if (length(available) == 0) stop("No requested features found in qc_flagged")

  long <- tidyr::pivot_longer(qc_flagged,
                               cols = dplyr::all_of(available),
                               names_to = "metric", values_to = "value")

  ggplot2::ggplot(long, ggplot2::aes(x = metric, y = value,
                                        color = flagged)) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
    ggplot2::scale_color_manual(values = c(`FALSE` = "grey40",
                                              `TRUE`  = "red")) +
    ggplot2::facet_wrap(~ metric, scales = "free", ncol = 3) +
    ggplot2::labs(x = NULL, y = "Value",
                   title = "Per-sample QC features",
                   subtitle = "Red = sample flagged by QC criteria") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom",
                    axis.text.x = ggplot2::element_blank())
}

# Small helper used above
`%||%` <- function(a, b) if (is.null(a)) b else a
