# Combine biomarker templates with case/control differential spec

Wraps templates in the group-level behavior spec used by
[`simulate_case_control_cohort()`](https://mljaniczek.github.io/classydms/reference/simulate_case_control_cohort.md).
Prevalences and intensity multipliers may be scalars (broadcast to every
biomarker) or vectors of length matching `nrow(templates)`.

## Usage

``` r
build_biomarker_spec(
  templates,
  case_prevalence,
  control_prevalence,
  case_intensity_mult = 1,
  control_intensity_mult = 1
)
```

## Arguments

- templates:

  tibble returned by
  [`holdout_biomarkers()`](https://mljaniczek.github.io/classydms/reference/holdout_biomarkers.md).

- case_prevalence, control_prevalence:

  Bernoulli probability that each biomarker appears in a case or control
  sample. Scalar or length-matched vector.

- case_intensity_mult, control_intensity_mult:

  Multiplier applied to the template `intensity` when the biomarker
  fires in a case or control sample. Scalar or length-matched vector.

## Value

tibble with the template columns plus `case_prevalence`,
`control_prevalence`, `case_intensity_mult`, `control_intensity_mult`.
Suitable for `biomarkers = ...` in
[`simulate_case_control_cohort()`](https://mljaniczek.github.io/classydms/reference/simulate_case_control_cohort.md).
