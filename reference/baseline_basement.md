# Dust-threshold a sample below a fixed intensity

Zeros out pixels of `sample$Z` whose intensity is below `basement_thr`.
Saves the pre-thresholded matrix into `sample$Z_pretrim` so peak-width
estimation can later use the unthresholded flanks.

## Usage

``` r
baseline_basement(sample, basement_thr)
```

## Arguments

- sample:

  Sample list.

- basement_thr:

  Threshold below which pixels are zeroed.

## Value

Updated sample.
