# Universal-backdrop biomarker preset for catalog simulation

Two-tier sample composition. Anchor layer: catalog compounds with
prevalence at or above `anchor_prevalence_threshold` (default 0.9) are
treated as universal-backdrop and appear in every simulated sample at
their catalog prevalence (typically ~1.0). Biomarker layer:
`n_biomarkers` non-anchor compounds are picked from the top intensity
quantile and given case/control Bernoulli prevalence overrides,
providing the classifier signal.

## Usage

``` r
catalog_biomarkers_universal_backdrop(
  catalog,
  n_biomarkers = 5L,
  anchor_prevalence_threshold = 0.9,
  biomarker_min_prevalence = 0.05,
  min_intensity_quantile = 0.95,
  case_prevalence = 0.9,
  control_prevalence = 0.1,
  case_intensity_mult = 1.2,
  control_intensity_mult = 1,
  seed = NULL
)
```

## Arguments

- catalog:

  A `peak_catalog` object.

- n_biomarkers:

  Number of biomarker compounds. Default `5L`.

- anchor_prevalence_threshold:

  Compounds with catalog prevalence `>=` this are anchors and excluded
  from the biomarker pool. Default `0.9`. On a universal catalog (built
  via
  [`build_universal_catalog()`](https://mljaniczek.github.io/classydms/reference/build_universal_catalog.md))
  this typically selects the ~150 compounds present across all sources.

- biomarker_min_prevalence:

  Minimum catalog prevalence for a non-anchor compound to be eligible as
  a biomarker. Default `0.05` — screens out ultra-rare compounds that
  may be jitter artifacts.

- min_intensity_quantile:

  Intensity quantile threshold applied to the compound's median observed
  intensity. Default `0.95`.

- case_prevalence, control_prevalence:

  Bernoulli firing probability of each biomarker in case / control
  samples. Defaults: `0.9` / `0.1`.

- case_intensity_mult, control_intensity_mult:

  Amplitude multiplier applied to biomarker peaks. Defaults: `1.2` /
  `1.0`.

- seed:

  Optional RNG seed.

## Value

A biomarker spec tibble (as
[`build_biomarker_spec_catalog()`](https://mljaniczek.github.io/classydms/reference/build_biomarker_spec_catalog.md))
plus a `"n_anchors"` attribute recording how many anchor compounds were
identified. Pass to
[`simulate_case_control_from_catalog()`](https://mljaniczek.github.io/classydms/reference/simulate_case_control_from_catalog.md).

## Details

Non-biomarker, non-anchor compounds keep their catalog empirical
prevalence (heterogeneous — some rare, some moderately common).

Biomarkers are picked from BELOW the anchor threshold so the signal
stands out against the invariant backdrop. Filtering to top
`min_intensity_quantile` (default 0.95) ensures they are strong enough
to be reliably detected.
