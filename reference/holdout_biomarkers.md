# Hold out K observed peaks as biomarker templates

Uniformly at random draws `n_biomarkers` peak indices from
`peak_params`'s raw pools, extracts their locations and morphology as
biomarker templates, and returns a modified `peak_params` with those
peaks *and* any peaks within a protected radius removed from the raw
pools. This guarantees that background samples drawn from the returned
`background_peak_params` will not include peaks that land near biomarker
locations, so injected biomarkers are truly additive signal rather than
contaminated by nearby background peaks.

## Usage

``` r
holdout_biomarkers(
  peak_params,
  n_biomarkers,
  H = NULL,
  W = NULL,
  min_sep_rt = 8L,
  min_sep_cv = 2L,
  location_jitter_rt = 2,
  location_jitter_cv = 1,
  jitter_safety_sd = 3,
  protect_radius_rt = NULL,
  protect_radius_cv = NULL,
  min_intensity_quantile = NULL,
  seed = NULL
)
```

## Arguments

- peak_params:

  Output of
  [`estimate_peak_params()`](https://mljaniczek.github.io/classydms/reference/estimate_peak_params.md);
  must contain `rt_loc_raw`, `cv_loc_raw`, `sigma_rt_raw`,
  `sigma_cv_raw`, `intensity_raw`.

- n_biomarkers:

  Integer; number of biomarker templates to draw.

- H, W:

  Optional integer image dimensions (pixels). When provided,
  `protect_radius_rt` and `protect_radius_cv` are automatically derived
  as `(min_sep + jitter_safety_sd * jitter) / dim`, which guarantees
  background peaks stay `min_sep + safety*jitter` pixels away from every
  biomarker location even after `location_jitter` is applied at
  simulation time. Strongly recommended over the fraction defaults,
  which were designed for a native ~1400 x 40 image and under-protect at
  smaller sizes.

- min_sep_rt, min_sep_cv:

  Minimum spatial separation between biomarker and background peak
  *centers* (in pixels). Defaults match the peak-detector defaults so
  the "protected" region is at least as large as the region a peak
  occupies.

- location_jitter_rt, location_jitter_cv:

  Per-peak location jitter SD (pixels) that will be applied at
  simulation time by
  [`simulate_case_control_cohort()`](https://mljaniczek.github.io/classydms/reference/simulate_case_control_cohort.md).
  Used when computing the derived `protect_radius_*` from `H`, `W`.

- jitter_safety_sd:

  Multiplier on `location_jitter_*` used in the derived protect radius
  formula. `3` covers ~99.7% of Normal draws so background peaks are
  essentially never jittered into a biomarker region.

- protect_radius_rt, protect_radius_cv:

  Radii (in fractions of image dimensions) around each biomarker inside
  which observed peaks are also removed from the background pool. When
  `H` (or `W`) is passed, these are computed automatically from the
  derived formula and any user-provided value is ignored. Backwards-
  compatible defaults (`15 / 1400` for RT, `3 / 40` for CV) are used
  only when neither `H` nor an explicit fraction is provided.

- min_intensity_quantile:

  Optional numeric in `[0, 1)`. When non-`NULL`, biomarker templates are
  drawn only from pool peaks whose `intensity_raw` is at or above the
  given quantile of the pool's intensity distribution — e.g. `0.90`
  restricts selection to the top-10% brightest peaks. Useful for
  sanity-check scenarios where you want the injected biomarker to be
  visually obvious even in a very sparse image. `NULL` (default) samples
  uniformly from the whole pool.

- seed:

  Optional RNG seed for reproducible template selection.

## Value

A list with:

- biomarker_templates:

  tibble with columns `rt_loc`, `cv_loc`, `sigma_rt`, `sigma_cv`,
  `intensity` — one row per held-out peak.

- background_peak_params:

  modified `peak_params` with the held-out peaks and any peaks within
  the protected radius stripped from every `*_raw` vector.
  `n_peaks_detected` is updated.

- n_pool_before:

  pool size before hold-out.

- n_pool_after:

  pool size after hold-out.
