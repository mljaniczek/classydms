# classydms: report rendering helpers

#' Render the peak-parameters diagnostic report
#'
#' Renders an HTML report summarizing the output of
#' [estimate_peak_params()], including per-sample peak counts, spatial
#' hotspot maps, sample-level clustering and shift analysis with
#' before-and-after global-alignment visualization, and per-class
#' distributions of peak widths, intensities, and background noise. When
#' class labels are supplied, class comparisons are reported alongside
#' Benjamini-Hochberg-adjusted p-values.
#'
#' `peak_params` must contain the `sample_idx_raw` field (parallel to
#' `rt_loc_raw`), which was added in classydms after commit
#' `Add sample_idx_raw and sample_names to estimate_peak_params`. Older
#' outputs of [estimate_peak_params()] will be rejected with an
#' informative error.
#'
#' Deferred to a follow-up (not covered here): k-means clustering on
#' per-sample shift vectors to detect unlabeled batch structure; anchor
#' peak analysis; sporadic-peak analysis; volcano plot; and the
#' downstream `build_peak_catalog()` function that this report is
#' designed to inform.
#'
#' @param peak_params Output of [estimate_peak_params()], containing the
#'   `sample_idx_raw` and `sample_names` fields.
#' @param y Optional named integer or factor vector of class labels
#'   (e.g. `2 = disease, 1 = control`). Names must match
#'   `peak_params$sample_names`. When `NULL`, class-comparison sections
#'   are skipped and the report is rendered label-agnostic.
#' @param batch Optional named vector (character or factor) of batch
#'   labels used to color the sample PCA. Names must match
#'   `peak_params$sample_names`. `NULL` skips batch coloring.
#' @param cohort_name Short label displayed in the report title
#'   (e.g. `"CHF"` or `"AKI"`).
#' @param output_file Output HTML filename.
#' @param output_dir Directory to render into. Defaults to the directory
#'   part of `output_file`; created if missing.
#' @param quiet Whether to suppress rmarkdown chatter (default `TRUE`).
#'
#' @return Path to the rendered HTML file (invisibly).
#' @export
render_peak_params_report <- function(peak_params,
                                       y = NULL,
                                       batch = NULL,
                                       cohort_name = "cohort",
                                       output_file = "peak_params_report.html",
                                       output_dir  = NULL,
                                       quiet = TRUE) {
  if (is.null(peak_params$sample_idx_raw)) {
    stop("peak_params$sample_idx_raw not found. Re-run estimate_peak_params ",
         "on a version of classydms that emits sample_idx_raw (added in the ",
         "'Add sample_idx_raw' commit).")
  }
  n_raw <- length(peak_params$rt_loc_raw)
  if (length(peak_params$sample_idx_raw) != n_raw) {
    stop("sample_idx_raw length (", length(peak_params$sample_idx_raw),
         ") does not match rt_loc_raw length (", n_raw, ").")
  }
  needed <- c("rmarkdown", "knitr", "ggplot2", "dplyr", "tibble")
  missing <- needed[!vapply(needed,
                            requireNamespace, logical(1),
                            quietly = TRUE)]
  if (length(missing)) {
    stop("The following packages are required to render the peak-params ",
         "report but are not installed: ",
         paste(missing, collapse = ", "), ".")
  }

  rmd_path <- system.file("reports", "peak_params_report.Rmd",
                          package = "classydms")
  if (!nzchar(rmd_path)) {
    # Fall back to devtools::load_all() layout when the package is not
    # installed but is being iterated on locally.
    dev_path <- file.path("inst", "reports", "peak_params_report.Rmd")
    if (file.exists(dev_path)) {
      rmd_path <- normalizePath(dev_path)
    } else {
      stop("peak_params_report.Rmd not found via system.file() or the ",
           "development inst/ path.")
    }
  }

  if (is.null(output_dir)) {
    output_dir <- dirname(normalizePath(output_file, mustWork = FALSE))
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  out_path <- rmarkdown::render(
    input       = rmd_path,
    output_file = basename(output_file),
    output_dir  = output_dir,
    params = list(
      peak_params = peak_params,
      y           = y,
      batch       = batch,
      cohort_name = cohort_name
    ),
    envir = new.env(parent = globalenv()),
    quiet = quiet
  )
  invisible(out_path)
}
