# Merge multiple peak catalogs into a single unified catalog

Combines compounds across catalogs by spatial proximity. Overlapping
compounds (within `(match_radius_rt, match_radius_cv)` of each other
across catalogs) are merged into a single entry: observations are
concatenated, and prevalence is recomputed as a weighted average by
cohort size. Compounds unique to a single catalog are added as new
entries.

## Usage

``` r
merge_catalogs(
  catalog_list,
  names = NULL,
  match_radius_rt = NULL,
  match_radius_cv = NULL
)
```

## Arguments

- catalog_list:

  List of `peak_catalog` objects.

- names:

  Optional character vector of length `length(catalog_list)` naming each
  source. Default `"catalog_1"`, `"catalog_2"`, etc.

- match_radius_rt, match_radius_cv:

  Cross-catalog matching radii. Defaults to the minimum `eps_rt` /
  `eps_cv` across the source catalogs — using a larger radius than the
  source eps would cause within-source compounds to fold into each other
  during the merge.

## Value

A `peak_catalog` object with the merged compounds.

## Details

The merged catalog has an extra `per_source_prevalence` list-column
recording per-input-catalog prevalence for each merged compound — useful
for tracking which sources contributed which compounds.
