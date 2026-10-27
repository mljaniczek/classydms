# Sanity-max scenario for catalog-based case/control simulation

Analogous to
[`biomarkers_sanity_max()`](https://mljaniczek.github.io/classydms/reference/biomarkers_sanity_max.md)
for the peak_params flow. Picks 5 catalog compounds whose intensity is
in the top 5% of compound median intensities and whose catalog
prevalence is at least `min_base_prevalence` (default 0.5). Sets case
prevalence to 1.0, control prevalence to 0.0, and case_intensity_mult to
10x. Designed to test that the full catalog-based pipeline can recover
an obviously injected signal.

## Usage

``` r
catalog_biomarkers_sanity_max(catalog, min_base_prevalence = 0.5, seed = NULL)
```

## Arguments

- catalog:

  A `peak_catalog` object.

- min_base_prevalence:

  Only pick biomarkers whose catalog prevalence is \>= this. Default
  `0.5`.

- seed:

  Optional RNG seed.

## Value

A biomarker spec tibble ready for
[`simulate_case_control_from_catalog()`](https://mljaniczek.github.io/classydms/reference/simulate_case_control_from_catalog.md).
