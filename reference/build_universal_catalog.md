# Build a universal catalog by pooling all samples from multiple cohorts

Unlike
[`merge_catalogs()`](https://mljaniczek.github.io/classydms/reference/merge_catalogs.md),
which reconciles already-built per-cohort catalogs, this function pools
raw peaks from all cohorts BEFORE clustering. Every peak from every
sample contributes directly to catalog construction, so systematic
per-cohort offsets get absorbed into the compound centroids rather than
producing duplicate entries.

## Usage

``` r
build_universal_catalog(
  pp_list,
  eps_rt = 0.01,
  eps_cv = 0.05,
  min_support_frac = 0.1,
  min_cluster_size = 2L,
  coord_mode = c("fraction", "physical"),
  align = FALSE
)
```

## Arguments

- pp_list:

  Named list of `peak_params` objects.

- eps_rt, eps_cv:

  Catalog clustering eps (fraction of image).

- min_support_frac:

  Support threshold. Applies to the combined sample set — a compound
  must span at least this fraction of the total N samples.

- min_cluster_size:

  Minimum peaks per cluster.

- align:

  Optional cross-cohort drift correction. `FALSE` (default) skips
  alignment. `TRUE` calls
  [`align_across_cohorts()`](https://mljaniczek.github.io/classydms/reference/align_across_cohorts.md)
  on `pp_list` with sensible defaults (`eps_rt_loose = 3 * eps_rt`,
  `eps_cv_loose = 3 * eps_cv`, `min_support_frac = min_support_frac`,
  `iterate = TRUE`, `max_iter = 3L`), then builds the universal catalog
  on the aligned cohorts. Corrects BOTH intra-cohort and inter-cohort
  drift in one step. Pass a named list to override any of those. Setting
  `align = TRUE` is a defensible answer to reviewer concerns about
  systematic per-sample or per-batch drift; the compactness-delta table
  in `compare_all_catalogs.Rmd` shows whether it helped on your cohorts.

## Value

A `peak_catalog` object with a `per_source_prevalence` list-column on
`compounds` giving cohort-level prevalence, and `source_names` /
`source_n_samples` fields in `parameters`.

## Details

Each catalog compound gets an extra list-column `per_source_prevalence`
recording the fraction of each input cohort's samples in which the
compound was observed. This mirrors
[`merge_catalogs()`](https://mljaniczek.github.io/classydms/reference/merge_catalogs.md)'s
output for downstream compatibility.
