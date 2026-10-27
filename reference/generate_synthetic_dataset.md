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
  size_jitter = 0.15,
  normalize = TRUE,
  dust_threshold = 0,
  location_mode = c("empirical", "marginal"),
  location_jitter_rt = 2,
  location_jitter_cv = 1,
  attribute_mode = c("joint", "marginal"),
  noise_scale = 1
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

  Per-peak random log-normal scale factor SD for size diversity. Default
  `0.15` (multiplicative variation of roughly x1.16 / :1.16), which
  matches typical GC-DMS peak-size variability and is the value used by
  the getting-started vignette. Larger values regularize more but
  produce peaks noticeably larger or smaller than any observed peak,
  which is not biologically plausible; use only if pretraining stability
  requires it.

- normalize:

  Whether to apply
  [`normalize_sample()`](https://mljaniczek.github.io/classydms/reference/normalize_sample.md)
  (default TRUE).

- dust_threshold:

  If \> 0, zero out pixels below this value in both the clean and noisy
  synthetic images BEFORE normalization. This matches the
  dust-thresholding step applied to real data via
  [`baseline_basement()`](https://mljaniczek.github.io/classydms/reference/baseline_basement.md)
  and produces synthetic samples whose sparsity profile matches real
  samples. Default 0 (no dust thresholding, backward-compatible). To
  match real-data preprocessing, set this to the same `basement_thr` you
  passed to
  [`baseline_basement()`](https://mljaniczek.github.io/classydms/reference/baseline_basement.md).

- location_mode:

  How to sample peak (RT, CV) locations. `"empirical"` (default,
  recommended) draws each peak's location by picking a random observed
  peak location from the real cohort (from `params$rt_loc_raw` /
  `params$cv_loc_raw`) and adding small per-peak jitter. This produces
  synthetic images whose peaks cluster at cohort-typical hotspots rather
  than scattering uniformly. `"marginal"` uses independent Normal draws
  on each axis with mean/SD from `params$rt_loc$mean`/`sd` — the
  previous behavior, kept for backward compatibility.

- location_jitter_rt, location_jitter_cv:

  SD of per-peak location jitter (in pixels) added on top of the
  empirical draw. Represents instrument reproducibility. Default 2
  pixels RT, 1 pixel CV.

- attribute_mode:

  How to sample per-peak morphology (sigma_rt, sigma_cv, intensity).
  `"joint"` (default, recommended) draws an index `j` into
  `params$rt_loc_raw` and reuses the observed
  `(sigma_rt_raw[j], sigma_cv_raw[j], intensity_raw[j])` triple. This
  preserves the empirical joint distribution: wider peaks are drawn with
  correspondingly higher intensity, matching the physical coupling in
  real GC-DMS data. When `location_mode == "empirical"` the same `j` is
  used for location and morphology; otherwise an independent index is
  drawn. `"marginal"` reverts to independent log-normal draws for each
  attribute (the pre-`attribute_mode` behavior; wider peaks and higher
  intensity are uncorrelated).

- noise_scale:

  Multiplier applied to `params$noise$sd` at generation time. Default
  `1.0` uses the estimated real-cohort noise level. Larger values
  corrupt the input more aggressively, which is a useful regularization
  knob for the denoising autoencoder (Vincent et al. 2008): too little
  corruption lets the encoder pass the input through with minimal
  compression; too much destroys the signal. A typical sweep is
  `noise_scale = c(1, 3, 5, 10)`, looking at real-validation MSE for the
  sweet spot.

## Value

List with `clean` and `noisy` 4-D torch tensors.
