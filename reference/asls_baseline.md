# Asymmetric Least Squares (AsLS) baseline correction

Iteratively estimates a smooth baseline that lies under the data using
Asymmetric Least Squares (AsLS) baseline correction

## Usage

``` r
asls_baseline(y, lambda = 1e+06, p = 0.01, niter = 10)
```

## Arguments

- y:

  Numeric vector (one RT trace).

- lambda:

  Smoothness penalty (default 1e6).

- p:

  Asymmetry parameter (default 0.01).

- niter:

  Number of iterations (default 10).

## Value

List with elements `z` (baseline) and `corrected` (`y - z`).

## Details

Iteratively estimates a smooth baseline that lies under the data using
asymmetric weights so peaks are not absorbed.
