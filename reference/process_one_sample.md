# Run preprocessing on one sample (list with `$Z`, optional `$time`, `$cv`)

Wraps
[`preprocess_matrix_rt()`](https://mljaniczek.github.io/classydms/reference/preprocess_matrix_rt.md)
and attaches the result back to the sample.

## Usage

``` r
process_one_sample(
  sample,
  sg_p = 3,
  sg_n = 21,
  asls_lambda = 1e+06,
  asls_p = 0.01,
  asls_niter = 10,
  clip_floor = 0
)
```

## Arguments

- sample:

  A list with at least an element `Z` (RT x CV matrix).

- sg_p, sg_n:

  Savitzky-Golay parameters.

- asls_lambda, asls_p, asls_niter:

  AsLS parameters.

- clip_floor:

  Optional floor (default 0).

## Value

The sample with `Z_raw`, `Z_smooth`, and updated `Z`.
