# Per-batch on-the-fly augmentation (RT/CV shifts + intensity scaling + noise)

Per-batch on-the-fly augmentation (RT/CV shifts + intensity scaling +
noise)

## Usage

``` r
augment_batch(
  x,
  rt_shift = 5,
  cv_shift = 2,
  intensity_scale = c(0.9, 1.1),
  noise_sd = 0,
  device = x$device
)
```
