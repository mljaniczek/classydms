# Pick a set of catalog compounds as biomarkers

Selects `n_biomarkers` compound_ids uniformly at random from the
catalog, optionally filtered by prevalence range or intensity quantile.
Analogous to
[`holdout_biomarkers()`](https://mljaniczek.github.io/classydms/reference/holdout_biomarkers.md)
for the peak_params flow.

## Usage

``` r
pick_catalog_biomarkers(
  catalog,
  n_biomarkers,
  min_prevalence = 0,
  max_prevalence = 1,
  min_intensity_quantile = NULL,
  seed = NULL
)
```

## Arguments

- catalog:

  A `peak_catalog` object.

- n_biomarkers:

  Number of biomarker compound_ids to pick.

- min_prevalence, max_prevalence:

  Only consider compounds whose catalog prevalence falls in this range.
  Defaults `[0, 1]` (any). For "rare-in-cohort" biomarkers set e.g.
  `max = 0.3`; for "common-in-cohort" set `min = 0.7`.

- min_intensity_quantile:

  Optional threshold on the median of each compound's `intensity_obs`.
  Same idea as in
  [`holdout_biomarkers()`](https://mljaniczek.github.io/classydms/reference/holdout_biomarkers.md).

- seed:

  Optional RNG seed.

## Value

Integer vector of length `n_biomarkers` — the selected compound_ids.

## Details

Unlike the peak_params version, no "protect radius" is needed — catalog
compounds already sit at discrete positions with no overlap between
them, so background compounds (all the non-selected entries) can't land
in biomarker regions.
