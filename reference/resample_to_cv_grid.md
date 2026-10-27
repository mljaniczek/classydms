# Resample one sample onto a CV grid

Interpolates each RT row of `$Z` from the sample's own `$cv` axis onto
`target_cv`. Columns outside the sample's CV span are filled with 0.
Needed because acquisitions can differ in CV sampling (e.g. 100 vs 101
columns over the same -40..15 V span).

## Usage

``` r
resample_to_cv_grid(sample, target_cv, method = c("linear", "constant"))
```

## Arguments

- sample:

  A sample list with `$cv` and `$Z` (`ncol(Z) == length(cv)`).

- target_cv:

  Increasing numeric vector of CV positions (volts).

- method:

  `"linear"` (default) or `"constant"` (nearest column).

## Value

The sample with `$Z` (`nrow(Z) x length(target_cv)`) and `$cv` set to
`target_cv`.
