# Simulate a case / control cohort with known injected biomarkers

Generates a synthetic cohort of `n_cases` case samples and `n_controls`
control samples. Each sample gets:

1.  A background drawn from `background_peak_params` via
    [`generate_one_synthetic()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic.md)
    with default `location_mode = "empirical"` and
    `attribute_mode = "joint"` (add_noise = FALSE at this stage; noise
    is added after biomarkers).

2.  For each biomarker row: a Bernoulli(`case_prevalence` or
    `control_prevalence`) draw decides whether the biomarker "fires" in
    this sample. If it fires, a 2D Gaussian peak with the biomarker's
    `(rt_loc, cv_loc, sigma_rt, sigma_cv)` and intensity
    `intensity * (case_intensity_mult | control_intensity_mult)` is
    added, plus per-peak jitter matching production
    (`location_jitter_rt/cv`, `size_jitter`).

3.  Folded-Normal background noise (from
    `background_peak_params$noise`).

## Usage

``` r
simulate_case_control_cohort(
  background_peak_params,
  biomarkers,
  n_cases,
  n_controls,
  H,
  W,
  size_jitter = 0.15,
  location_jitter_rt = 2,
  location_jitter_cv = 1,
  add_noise = TRUE,
  seed = NULL
)
```

## Arguments

- background_peak_params:

  Modified `peak_params` returned by
  [`holdout_biomarkers()`](https://mljaniczek.github.io/classydms/reference/holdout_biomarkers.md)
  (or the same, minus the biomarker peaks).

- biomarkers:

  Output of
  [`build_biomarker_spec()`](https://mljaniczek.github.io/classydms/reference/build_biomarker_spec.md)
  (or a compatible tibble with the same columns).

- n_cases, n_controls:

  Number of samples in each group.

- H, W:

  Output image dimensions.

- size_jitter, location_jitter_rt, location_jitter_cv:

  Passed to the background generator; also applied to biomarker peaks so
  their morphology varies per sample.

- add_noise:

  Whether to add background noise (default TRUE).

- seed:

  Optional RNG seed.

## Value

A list with:

- samples:

  named list of length `n_cases + n_controls`. Each element is
  `list(path, time, cv, Z)` matching the shape produced by
  [`load_dms_directory()`](https://mljaniczek.github.io/classydms/reference/load_dms_directory.md),
  so it can flow through the rest of the classydms pipeline unchanged.

- y:

  named integer vector: 2 = case, 1 = control.

- ground_truth:

  the biomarkers tibble.

- per_sample_firing:

  logical matrix (`n_samples` x `n_biomarkers`) recording which
  biomarkers actually fired in each sample, for post-hoc analysis.
