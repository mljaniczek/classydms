# Generate N synthetic pairs as torch tensors

Calls
[`generate_one_synthetic()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic.md)
N times, applies per-sample log-quantile normalization, and stacks the
result as torch tensors of shape `(N, 1, H, W)`.

## Usage

``` r
generate_synthetic_dataset(
  params,
  N,
  H,
  W,
  add_noise = TRUE,
  size_jitter = 0.6,
  normalize = TRUE
)
```

## Arguments

- params:

  Output of
  [`estimate_peak_params()`](https://mljaniczek.github.io/classydms/reference/estimate_peak_params.md).

- N:

  Number of synthetic samples to generate.

- H, W:

  Output image dimensions.

- add_noise:

  If TRUE, add spatially-varying background noise.

- size_jitter:

  Per-peak random log-normal scale factor SD for size diversity (default
  0.6).

- normalize:

  Whether to apply
  [`normalize_sample()`](https://mljaniczek.github.io/classydms/reference/normalize_sample.md)
  (default TRUE).

## Value

List with `clean` and `noisy` 4-D torch tensors.
