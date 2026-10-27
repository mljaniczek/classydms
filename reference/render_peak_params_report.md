# Render the peak-parameters diagnostic report

Renders an HTML report summarizing the output of
[`estimate_peak_params()`](https://mljaniczek.github.io/classydms/reference/estimate_peak_params.md),
including per-sample peak counts, spatial hotspot maps, sample-level
clustering and shift analysis with before-and-after global-alignment
visualization, and per-class distributions of peak widths, intensities,
and background noise. When class labels are supplied, class comparisons
are reported alongside Benjamini-Hochberg-adjusted p-values.

## Usage

``` r
render_peak_params_report(
  peak_params,
  y = NULL,
  batch = NULL,
  class_labels = NULL,
  cohort_name = "cohort",
  output_file = "peak_params_report.html",
  output_dir = NULL,
  quiet = TRUE
)
```

## Arguments

- peak_params:

  Output of
  [`estimate_peak_params()`](https://mljaniczek.github.io/classydms/reference/estimate_peak_params.md),
  containing the `sample_idx_raw` and `sample_names` fields.

- y:

  Optional named integer or factor vector of class labels (e.g.
  `2 = disease, 1 = control`). Names must match
  `peak_params$sample_names`. When `NULL`, class-comparison sections are
  skipped and the report is rendered label-agnostic.

- batch:

  Optional named vector (character or factor) of batch labels used to
  color the sample PCA. Names must match `peak_params$sample_names`.
  `NULL` skips batch coloring.

- class_labels:

  Optional named character vector remapping the values of `y` to display
  labels (e.g. `c("1" = "control", "2" = "disease")`). If provided,
  class comparisons and plot legends use the display labels; values in
  `y` not present in `names(class_labels)` are shown unchanged. Ignored
  when `y` is already a factor with meaningful character levels.

- cohort_name:

  Short label displayed in the report title (e.g. `"CHF"` or `"AKI"`).

- output_file:

  Output HTML filename.

- output_dir:

  Directory to render into. Defaults to the directory part of
  `output_file`; created if missing.

- quiet:

  Whether to suppress rmarkdown chatter (default `TRUE`).

## Value

Path to the rendered HTML file (invisibly).

## Details

`peak_params` must contain the `sample_idx_raw` field (parallel to
`rt_loc_raw`), which was added in classydms after commit
`Add sample_idx_raw and sample_names to estimate_peak_params`. Older
outputs of
[`estimate_peak_params()`](https://mljaniczek.github.io/classydms/reference/estimate_peak_params.md)
will be rejected with an informative error.

Deferred to a follow-up (not covered here): k-means clustering on
per-sample shift vectors to detect unlabeled batch structure; anchor
peak analysis; sporadic-peak analysis; volcano plot; and the downstream
[`build_peak_catalog()`](https://mljaniczek.github.io/classydms/reference/build_peak_catalog.md)
function that this report is designed to inform.
