# Level 1 (sanity check) biomarker scenario

1 biomarker drawn from the *top 10% brightest* pool peaks
(`min_intensity_quantile = 0.90`), extreme prevalence gap (95% cases vs
5% controls), 5x intensity multiplier in cases. Designed to be visually
obvious in a spot-check plot; case/control AUC should be well above 0.90
if the pipeline works. If Level 1 fails, something in the pipeline is
broken.

## Usage

``` r
biomarkers_level1(peak_params, H = NULL, W = NULL, seed = NULL)
```

## Arguments

- peak_params:

  Output of
  [`estimate_peak_params()`](https://mljaniczek.github.io/classydms/reference/estimate_peak_params.md).

- H, W:

  Optional image dimensions in pixels. When provided, forwarded to
  [`holdout_biomarkers()`](https://mljaniczek.github.io/classydms/reference/holdout_biomarkers.md)
  so the protect radius scales safely for the target image size.
  Strongly recommended — the default fraction-based protect radius
  under-protects on smaller images than ~1400 x 40 and can cause
  background peaks to leak into biomarker regions after location jitter.

- seed:

  Optional RNG seed for reproducible biomarker selection.

## Value

list(biomarkers, background_peak_params) — pass both to
[`simulate_case_control_cohort()`](https://mljaniczek.github.io/classydms/reference/simulate_case_control_cohort.md).
