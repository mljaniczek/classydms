# Quantitative comparison of synthetic vs real GC-DMS samples

Generates a batch of synthetic samples from `peak_params` and compares
their distributional summary statistics to those of real samples. Use
this to verify the synthetic generator's output resembles real spectra
before committing to a long pretraining run.

## Usage

``` r
synthetic_quality_check(
  Z_real_list,
  peak_params,
  H = NULL,
  W = NULL,
  n_synthetic = 50L,
  size_jitter = 0.15,
  normalize = FALSE,
  top_k = 10000L,
  dust_threshold = 0,
  location_mode = c("empirical", "marginal"),
  location_jitter_rt = 2,
  location_jitter_cv = 1,
  attribute_mode = c("joint", "marginal"),
  noise_scale = 1,
  seed = NULL
)
```

## Arguments

- Z_real_list:

  List of real preprocessed Z matrices. Typically the same `Z_raw` you
  would pass to
  [`estimate_peak_params()`](https://mljaniczek.github.io/classydms/reference/estimate_peak_params.md).

- peak_params:

  Output of
  [`estimate_peak_params()`](https://mljaniczek.github.io/classydms/reference/estimate_peak_params.md).

- H, W:

  Dimensions for synthetic generation. Default uses the dimensions of
  the first real sample.

- n_synthetic:

  Number of synthetic samples to generate.

- size_jitter:

  Same as in
  [`generate_one_synthetic()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic.md).

- normalize:

  If TRUE, apply
  [`normalize_sample()`](https://mljaniczek.github.io/classydms/reference/normalize_sample.md)
  to both real and synthetic matrices before comparing. Default FALSE
  compares on the raw (post-baseline-correction, pre-normalization)
  scale.

- top_k:

  Maximum peaks counted per sample for the n_peaks comparison via
  [`pick_peak_centers()`](https://mljaniczek.github.io/classydms/reference/pick_peak_centers.md).
  Set high enough to NOT cap detected peak counts (e.g., 10000) so the
  comparison reflects actual peak density rather than the cap. Default
  10000.

- dust_threshold:

  If \> 0, apply the same dust threshold to synthetic samples that was
  applied to real data. Pass the same value used in
  [`baseline_basement()`](https://mljaniczek.github.io/classydms/reference/baseline_basement.md)
  so synthetic and real are compared on the same preprocessed sparsity
  profile.

- location_mode:

  Passed through to
  [`generate_one_synthetic()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic.md).
  Default `"empirical"` samples peak locations from the real cohort's
  observed hotspots.

- location_jitter_rt, location_jitter_cv:

  Passed through to
  [`generate_one_synthetic()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic.md).

- attribute_mode:

  Passed through to
  [`generate_one_synthetic()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic.md).
  Default `"joint"` couples per-peak sigma_rt / sigma_cv / intensity by
  resampling the observed triple.

- seed:

  Optional RNG seed for reproducibility.

## Value

A list with:

- summary:

  data.frame of comparison metrics, one row per statistic.

- intensity_real:

  pooled positive pixel intensities from real samples (numeric vector).

- intensity_synth:

  same for synthetic samples.

- per_sample:

  data.frame of per-sample stats for both groups, useful for plotting.

## Details

Five comparisons are reported, each with a Kolmogorov-Smirnov D
statistic and p-value (small p-value indicates the real and synthetic
distributions differ):

1.  **pixel_intensity_pooled**: distribution of all positive pixel
    values, pooled across samples. Tests whether the intensity dynamic
    range matches.

2.  **total_signal**: per-sample sum of all pixels. Tests overall signal
    magnitude.

3.  **mean_positive**: per-sample mean of positive pixels. Tests typical
    peak amplitude scale.

4.  **sparsity**: per-sample fraction of zero pixels. Tests whether the
    proportion of background-vs-signal region matches.

5.  **n_peaks**: number of detected peaks per sample (via
    [`pick_peak_centers()`](https://mljaniczek.github.io/classydms/reference/pick_peak_centers.md)).
    Tests whether the peak density matches.

## See also

[`plot_synthetic_quality()`](https://mljaniczek.github.io/classydms/reference/plot_synthetic_quality.md)
for a visual companion.
