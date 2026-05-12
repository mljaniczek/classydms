# Compute a per-sample relative "dust" threshold

Returns 1% of the q-th quantile of positive values, used to identify
near-zero pixels that are likely numerical noise rather than real
signal. Returns 0 for samples with fewer than 50 positive values
(treated as junk).

## Usage

``` r
robust_eps(Z, q = 0.95, frac = 0.01)
```

## Arguments

- Z:

  A numeric matrix.

- q:

  Quantile to base the threshold on (default 0.95).

- frac:

  Fraction of the quantile to use (default 0.01).

## Value

A scalar threshold value.
