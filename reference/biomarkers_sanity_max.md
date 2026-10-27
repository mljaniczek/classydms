# Maximally-favorable ("sanity-max") biomarker scenario

Combines every lever we have in the direction of "make signal impossible
to miss":

- 5 biomarkers (collective signal, not one lone peak)

- Drawn from top 5% of pool by intensity

- Perfect presence / absence: `case_prevalence = 1.0`,
  `control_prevalence = 0.0`

- Case intensity multiplier = 10x

- Sparse background: `n_peaks$mean = 30` and `n_peaks$sd = 5`
  overwritten on the returned peak_params so synthetic samples have only
  ~30 background peaks rather than ~1300, making injected biomarkers
  dominate total image intensity

## Usage

``` r
biomarkers_sanity_max(peak_params, H = NULL, W = NULL, seed = NULL)
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

## Details

Intended as a **pipeline validation** scenario, not a difficulty
benchmark: if case/control AUC on this scenario is \< 0.95, something in
the pipeline is broken and Levels 1/2/3 are not worth running.
