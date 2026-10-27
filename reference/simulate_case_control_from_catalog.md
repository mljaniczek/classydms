# Simulate a case / control cohort from a peak catalog

For each sample: determine group (case or control), then generate via
[`generate_one_synthetic_from_catalog()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic_from_catalog.md)
with the biomarker spec's per-group prevalence and intensity overrides
applied. Non-biomarker compounds use the catalog's own prevalence in
both groups (background is shared). This directly matches the real-data
property that the *same* compound positions appear across most samples,
only differing in whether individual compounds fired.

## Usage

``` r
simulate_case_control_from_catalog(
  catalog,
  biomarkers,
  n_cases,
  n_controls,
  H,
  W,
  size_jitter = 0.15,
  location_jitter_rt = 2,
  location_jitter_cv = 1,
  add_noise = TRUE,
  noise_scale = 1,
  cofire = NULL,
  off_catalog_peaks = NULL,
  peak_shape = c("gaussian", "emg", "mix"),
  emg_tau_range = c(0.2, 0.8),
  emg_mix_prob = 0.3,
  seed = NULL
)
```

## Arguments

- catalog:

  A `peak_catalog` object.

- biomarkers:

  A tibble from
  [`build_biomarker_spec_catalog()`](https://mljaniczek.github.io/classydms/reference/build_biomarker_spec_catalog.md).

- n_cases, n_controls:

  Number of samples per group.

- H, W:

  Output image dimensions in pixels.

- size_jitter, location_jitter_rt, location_jitter_cv:

  Per-peak jitter parameters, forwarded to the underlying generator.

- add_noise:

  Whether to add background noise.

- noise_scale:

  Multiplier on catalog noise SD.

- cofire:

  Optional list of dependency groups modelling co-firing structure
  between compounds. Each element is a list with EITHER `anchor` (a
  single `compound_id`), `dependents` (integer vector of
  `compound_id`s), and `p_dependent_given_anchor` (scalar in `[0,1]` —
  probability each dependent fires when the anchor does), OR `members`
  (integer vector of `compound_id`s) and `p_joint` (scalar in `[0,1]` —
  probability the whole clique fires together). Applied as
  post-processing on top of the independent Bernoulli firing:

  - Dependency group: if anchor did not fire, every dependent that
    independently fired is REMOVED (its Gaussian subtracted). If the
    anchor fired, every dependent that did not fire independently gets
    an additional Bernoulli(`p_dependent_given_anchor`) roll; on success
    the compound fires and its Gaussian is added.

  - Symmetric clique: one joint Bernoulli(`p_joint`) roll decides the
    whole group. On success every non-fired member is added; on failure
    every fired member is removed. The final `per_sample_firing` matrix
    reflects the post-cofire state.

- off_catalog_peaks:

  Optional list controlling injection of random peaks NOT present in the
  catalog. Fields:

  - `n_cases`, `n_controls`: integer peaks to inject per case / control
    sample.

  - `rt_range`, `cv_range`: length-2 numerics in 0,1 giving the
    placement region as fractions of the padded `H x W` canvas (position
    = fraction \* H). On a centre-padded canvas the data block occupies
    only part of this range; see
    [`catalog_pixel_positions()`](https://mljaniczek.github.io/classydms/reference/catalog_pixel_positions.md).

  - `intensity_range`: length-2 numeric giving uniform bounds for the
    injected peak amplitude.

  - `exclusion_eps_rt`, `exclusion_eps_cv`: fractional radius; a
    candidate position is rejected if within this box of any catalog
    compound. Defaults `(0.02, 0.05)`. Sigmas are drawn from the pooled
    catalog sigma distribution so peak shapes match the catalog's
    morphology. Used for stress-testing an encoder / classifier on peaks
    the pretraining never saw. Peak positions are recorded on the
    returned `off_catalog_positions` list (one entry per sample).

- seed:

  Optional RNG seed.

## Value

List with `samples`, `y`, `ground_truth` (the biomarker spec),
`per_sample_firing` (logical matrix of `n_samples` x `n_biomarkers`;
whether each biomarker fired in each sample AFTER any cofire
post-processing), `cofire` (the cofire spec that was applied, or
`NULL`), and `off_catalog_positions` (list of length `n_samples`
recording injected positions, or `NULL`). `samples` matches
[`load_dms_directory()`](https://mljaniczek.github.io/classydms/reference/load_dms_directory.md)'s
shape so downstream pipeline steps work unchanged.

## Details

Compared to
[`simulate_case_control_cohort()`](https://mljaniczek.github.io/classydms/reference/simulate_case_control_cohort.md)
(the peak_params flow):

- Shared background: every sample draws from the same catalog compound
  positions rather than random-per-sample peak positions.

- No protect_radius: catalog compounds are already at discrete positions
  with no overlap between them.

- Biomarker specification is by `compound_id`, not by random draw from a
  pool.
