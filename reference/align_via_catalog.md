# Align samples via cluster-membership offsets

Alternative alignment approach that avoids the distance-matching
fragility of `align_peak_params()`. Steps:

## Usage

``` r
align_via_catalog(
  peak_params,
  eps_rt_loose = 0.03,
  eps_cv_loose = 0.08,
  min_support_frac = 0.1,
  min_samples_per_anchor = NULL,
  min_anchors_per_sample = 5L,
  iterate = TRUE,
  max_iter = 3L,
  verbose = TRUE
)
```

## Arguments

- peak_params:

  `peak_params` object with `sample_idx_raw`.

- eps_rt_loose, eps_cv_loose:

  Preliminary-catalog clustering eps. Must be **larger than the expected
  drift magnitude** so drifted samples' peaks join the same cluster as
  un-drifted samples' peaks. Defaults `(0.03, 0.08)` accommodate up to
  ~3% RT drift.

- min_support_frac:

  Support fraction for the preliminary catalog.

- min_samples_per_anchor:

  A compound must have at least this many unique samples to serve as an
  alignment anchor. `NULL` (default) sets it to
  `max(3, floor(0.30 * n_samples))` so only well-supported compounds
  contribute offsets.

- min_anchors_per_sample:

  A sample won't be shifted if fewer anchor compounds contain a peak
  from it.

- iterate:

  Whether to rebuild catalog on aligned peaks and re-align. Default
  `TRUE`.

- max_iter:

  Maximum iterations.

- verbose:

  Whether to print progress.

## Value

Updated `peak_params` with aligned `rt_loc_raw`, `cv_loc_raw`, and an
`alignment` field with method = `"cluster_membership"`, `shifts_rt`,
`shifts_cv`, `n_anchors_matched` (per sample), `n_anchors_total`,
`iterations`, `converged`.

## Details

1.  Build a preliminary catalog with **loose** eps (large enough to
    capture drifted peaks in the same cluster as un-drifted ones)

2.  For each anchor compound with sufficient sample support, each member
    peak has a known offset from the cluster centroid

3.  For each sample, aggregate its member-peak offsets across all anchor
    compounds — the median gives that sample's drift estimate

4.  Apply per-sample shifts to all peaks

5.  Optionally iterate: rebuild catalog on aligned peaks and refine

Why this is more robust than distance-matching: peak-to-anchor
correspondences are established by the preliminary clustering, not by
nearest-neighbor search. In dense peak clouds where many unrelated peaks
fall near each anchor, distance-matching gets confused;
cluster-membership does not.

Signals that this alignment is helping (after building the final catalog
on aligned peaks vs original): `n_compounds` decreases (fragments
consolidate), `n_singletons` decreases, `median_prevalence` increases,
and `n_anchor_90pct` increases.
