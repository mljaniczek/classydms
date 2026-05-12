# Savitzky-Golay smoothing + AsLS baseline correction (per CV column)

Smooths each CV column along RT with a Savitzky-Golay filter, then
subtracts an AsLS-estimated baseline.

## Usage

``` r
preprocess_matrix_rt(
  Z,
  sg_p = 3,
  sg_n = 21,
  asls_lambda = 1e+06,
  asls_p = 0.01,
  asls_niter = 10,
  clip_floor = 0
)
```

## Arguments

- Z:

  RT x CV intensity matrix.

- sg_p:

  Polynomial order for Savitzky-Golay (default 3).

- sg_n:

  Window length (default 21; rounded up to odd, capped at nrow).

- asls_lambda, asls_p, asls_niter:

  AsLS parameters.

- clip_floor:

  Optional floor for corrected values (default 0).

## Value

List with `Z_smooth` (smoothed only) and `Z_corrected` (smoothed +
baseline-corrected).
