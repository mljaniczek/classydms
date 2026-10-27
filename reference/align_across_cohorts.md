# Align samples across multiple cohorts jointly

Extends
[`align_via_catalog()`](https://mljaniczek.github.io/classydms/reference/align_via_catalog.md)
across cohorts. Merges all `peak_params` objects into one super-cohort
with globally-unique sample indices, aligns using cluster-membership
offsets, then splits the aligned peaks back into per-cohort
`peak_params`.

## Usage

``` r
align_across_cohorts(
  pp_list,
  eps_rt_loose = 0.04,
  eps_cv_loose = 0.1,
  min_support_frac = 0.1,
  min_samples_per_anchor = NULL,
  min_anchors_per_sample = 5L,
  iterate = TRUE,
  max_iter = 3L,
  verbose = TRUE
)
```

## Arguments

- pp_list:

  Named list of `peak_params` objects (one per cohort). All must have
  `sample_idx_raw`. Names become cohort labels in the returned alignment
  metadata.

- eps_rt_loose, eps_cv_loose:

  Preliminary catalog eps. Must be larger than the expected drift
  magnitude across cohorts. Default `(0.04, 0.10)` is looser than the
  intra-cohort version because inter-cohort drift is often larger.

- min_support_frac:

  Support threshold for the merged catalog. Kept low (default `0.10`)
  because merged sample counts are large.

- min_samples_per_anchor:

  Minimum unique samples an anchor compound must span. `NULL` (default)
  sets to `max(5, floor(0.20 * total_n_samples))`.

- min_anchors_per_sample:

  Minimum anchor compounds a sample must contribute to for a shift to be
  computed.

- iterate, max_iter, verbose:

  Same as
  [`align_via_catalog()`](https://mljaniczek.github.io/classydms/reference/align_via_catalog.md).

## Value

A named list of aligned `peak_params` objects (same names as `pp_list`).
Each carries an `alignment` field with method = `"cross_cohort"` plus
the per-sample shifts, iteration count, and the number of anchor
compounds used.

## Details

This corrects **cross-cohort drift** — the situation where samples
within any single cohort are well-aligned relative to each other, but
different cohorts' compound positions sit at systematically different
`(rt_loc, cv_loc)` because of instrument drift over months or between
studies. Intra-cohort alignment cannot detect or fix this because
there's no drift signal within each cohort's own samples.

When alignment is run across cohorts jointly, the preliminary catalog
sees compound positions across ALL sources at once, so the cluster
centroids reflect the average position across cohorts. Each cohort's
samples then get shifted to that shared reference, collapsing the
cross-cohort duplication in the universal catalog.
