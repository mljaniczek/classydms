# Build a case/control biomarker spec for catalog-based simulation

Wraps selected compound_ids with the per-group prevalence and intensity
overrides used by
[`simulate_case_control_from_catalog()`](https://mljaniczek.github.io/classydms/reference/simulate_case_control_from_catalog.md).

## Usage

``` r
build_biomarker_spec_catalog(
  compound_ids,
  case_prevalence,
  control_prevalence,
  case_intensity_mult = 1,
  control_intensity_mult = 1
)
```

## Arguments

- compound_ids:

  Integer vector of compound_ids designated as biomarkers (see
  [`pick_catalog_biomarkers()`](https://mljaniczek.github.io/classydms/reference/pick_catalog_biomarkers.md)).

- case_prevalence, control_prevalence:

  Bernoulli probability of the biomarker firing in a case / control
  sample. Scalar (broadcast) or length-matched vector.

- case_intensity_mult, control_intensity_mult:

  Amplitude multiplier applied to the compound's bootstrapped intensity
  in case / control samples. Scalar or length-matched vector.

## Value

A tibble with columns `compound_id`, `case_prevalence`,
`control_prevalence`, `case_intensity_mult`, `control_intensity_mult`.
Pass to
[`simulate_case_control_from_catalog()`](https://mljaniczek.github.io/classydms/reference/simulate_case_control_from_catalog.md).
