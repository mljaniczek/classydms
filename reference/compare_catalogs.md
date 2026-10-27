# Compare two peak catalogs by compound proximity

Matches each compound in `catalog_a` to its nearest neighbor in
`catalog_b` by spatial proximity in fraction-of-image `(rt_loc, cv_loc)`
space. Compounds are "matched" if the nearest neighbor is within
`(match_radius_rt, match_radius_cv)`. Returns three data frames: shared
compounds (matched to each other), compounds only in A, and compounds
only in B. Also reports the per-catalog prevalence for shared compounds
so you can see whether the same compound expresses at different rates
between cohorts.

## Usage

``` r
compare_catalogs(
  catalog_a,
  catalog_b,
  name_a = "a",
  name_b = "b",
  match_radius_rt = 0.01,
  match_radius_cv = 0.05
)
```

## Arguments

- catalog_a, catalog_b:

  `peak_catalog` objects.

- name_a, name_b:

  Short labels used as column-name suffixes in the `shared` output.
  Default `"a"` and `"b"`.

- match_radius_rt, match_radius_cv:

  Maximum distance to count a compound as matched. Defaults
  `(0.01, 0.05)` — the same defaults as `build_peak_catalog`.

## Value

List with:

- shared:

  tibble of matched compound pairs with `compound_id_<a>`,
  `compound_id_<b>`, `rt_loc_<a>`, `rt_loc_<b>`, `cv_loc_<a>`,
  `cv_loc_<b>`, `prevalence_<a>`, `prevalence_<b>`, `prevalence_delta`
  (b − a), and `spatial_distance` (Euclidean in the fraction space).

- only_a:

  tibble of compounds present in A but not matched in B.

- only_b:

  tibble of compounds present in B but not matched in A.

- summary:

  one-row tibble of counts + fractions.
