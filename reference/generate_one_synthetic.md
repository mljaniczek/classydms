# Generate one synthetic GC-DMS image

Generates a single synthetic image as a sum of 2-D Gaussian peaks
(sampled from the empirical peak distributions in `params`) plus
optional background noise.

## Usage

``` r
generate_one_synthetic(params, H, W, add_noise = TRUE, size_jitter = 0.6)
```

## Arguments

- params:

  Output of
  [`estimate_peak_params()`](https://mljaniczek.github.io/classydms/reference/estimate_peak_params.md).

- H, W:

  Output image dimensions.

- add_noise:

  If TRUE, add spatially-varying background noise.

- size_jitter:

  Per-peak random log-normal scale factor SD for size diversity (default
  0.6).

## Value

List with `clean` and `noisy` matrices (`clean == noisy` if
`add_noise = FALSE`).
