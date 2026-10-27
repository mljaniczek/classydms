# Merge multiple peak_params objects into a single super-cohort

Concatenates all `_raw` vectors across the input peak_params objects
with globally-unique `sample_idx_raw`. Useful for building a universal
catalog from all cohorts jointly.

## Usage

``` r
merge_peak_params(pp_list)
```

## Arguments

- pp_list:

  Named list of `peak_params` objects. All must have `sample_idx_raw`.
  Names become cohort labels.

## Value

A single `peak_params`-shaped list with concatenated `_raw` vectors,
renumbered `sample_idx_raw`, `sample_cohort` giving per-sample cohort
membership, and `cohort_names` (input names).

## Details

Each returned merged sample's original cohort is recorded in the
`sample_cohort` field, which is a character vector of length `n_samples`
giving the input list's name for each sample.
