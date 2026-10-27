# Generate one three-layer synthetic sample for encoder pretraining

Builds a synthetic GC-DMS image by stacking three layers of peaks:

## Usage

``` r
generate_noisy_pretrain_sample(
  catalog,
  H,
  W,
  anchor_ids,
  variable_config = NULL,
  contamination_config = NULL,
  add_noise = TRUE,
  size_jitter = 0.15,
  location_jitter_rt = 2,
  location_jitter_cv = 1,
  noise_scale = 1,
  anchor_dropout = 0,
  anchor_jitter_variability = NULL,
  variable_source_peak_params = NULL,
  peak_shape = c("gaussian", "emg", "mix"),
  emg_tau_range = c(0.2, 0.8),
  emg_mix_prob = 0.3,
  anchor_draw = c("all", "prevalence")
)
```

## Arguments

- catalog:

  A `peak_catalog` object.

- H, W:

  Output image dimensions in pixels.

- anchor_ids:

  Integer vector of compound_ids to force-fire. Typically the ~100
  anchors from the universal catalog with prevalence \>= 0.9.

- variable_config:

  Optional list controlling the variable layer. Fields:

  - `draw`: how compounds are chosen from the eligible pool. `"uniform"`
    (default, the original behaviour) draws `n_per_sample` of them
    uniformly without replacement, so every eligible compound appears in
    about `n / length(eligible)` of images whatever its own prevalence.
    `"prevalence"` fires each eligible compound independently with
    probability `min(1, prevalence * prevalence_scale)` (the same
    Bernoulli rule as
    [`generate_one_synthetic_from_catalog()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic_from_catalog.md)),
    so the expected count is the sum of those probabilities and
    per-compound frequencies match the catalog. Needs a `prevalence`
    column.

  - `prevalence_scale`: non-negative multiplier on catalog prevalence in
    `"prevalence"` mode. Default `1`. Ignored in `"uniform"` mode.

  - `n_per_sample`: scalar (fixed count) or length-2 numeric (uniform
    range). In `"uniform"` mode, how many variable compounds to draw;
    default `c(30, 80)`. In `"prevalence"` mode, an optional per-sample
    cap: when more compounds fire than the cap, a uniform subset of that
    size is kept. Default `NULL` (no cap); a cap below the expected
    count biases frequencies back toward uniform.

  - `prevalence_range`: length-2 numeric giving the catalog prevalence
    range eligible for the variable pool. Default `c(0.05, 0.5)`.

  - `intensity_scale`: length-2 numeric giving amplitude multiplier
    bounds. Default `c(0.3, 1.0)` (variable compounds slightly dimmer
    than anchors on average).

  - `exclude_anchors`: `TRUE` (default) — remove every compound in
    `anchor_ids` from the variable pool so anchors don't get drawn
    twice. This uses the full `anchor_ids`, including anchors that
    `anchor_draw = "prevalence"` or `anchor_dropout` silenced this
    sample, so a dropped anchor is not re-drawn as a variable compound.
    Pass `NULL` to skip this layer.

- contamination_config:

  Optional list controlling contamination injection. Same fields as
  `off_catalog_peaks` in
  [`simulate_case_control_from_catalog()`](https://mljaniczek.github.io/classydms/reference/simulate_case_control_from_catalog.md)
  except with `n_per_sample` instead of `n_cases` / `n_controls`.
  Default field values: `n_per_sample = c(5, 30)`, `rt_range` /
  `cv_range` = 5-95% of the catalog's (median) unpadded data block
  expressed as canvas fractions (5-95% of the canvas for catalogs
  without data dims), `intensity_range = c(0.1, 0.5)`,
  `exclusion_eps_rt = 0.02`, `exclusion_eps_cv = 0.05`. Pass `NULL` to
  skip this layer. The exclusion boxes surround every catalog compound
  position, not only the compounds placed in this sample, so they do not
  depend on the anchor / variable draw. Note that these defaults are
  much brighter and more exclusive than real off-catalog peaks: on a
  606-sample universal catalog, real off-catalog peaks had median
  normalised intensity about 0.014 (90th percentile about 0.042) against
  the default `intensity_range = c(0.1, 0.5)`, and the default exclusion
  boxes confine contamination to the gaps between compounds (only 61 of
  96 CV columns reachable on a 1504 x 96 grid). For realistic clutter
  consider `intensity_range` near `c(0.005, 0.05)` and smaller boxes
  (e.g. `exclusion_eps_rt = 0.002`, `exclusion_eps_cv = 0.01`).

- add_noise:

  Whether to add folded-Normal background noise on top of the three
  layers.

- size_jitter, location_jitter_rt, location_jitter_cv:

  Per-peak jitter parameters (same convention as
  [`generate_one_synthetic_from_catalog()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic_from_catalog.md)).

- noise_scale:

  Multiplier on catalog noise SD.

- anchor_dropout:

  Per-sample Bernoulli probability of dropping each anchor before
  rendering. `0` (default) fires all anchors.

- anchor_jitter_variability:

  Optional length-2 range for a per-sample scalar multiplier on
  `location_jitter_rt/cv`. `NULL` (default) uses the fixed jitter values
  as-is.

- variable_source_peak_params:

  Optional `peak_params` bundle. When provided, the variable layer is
  rendered via
  [`generate_one_synthetic()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic.md)
  on this bundle (real peaks, empirical locations, joint attributes)
  instead of the catalog-clustered render. Mutually exclusive with
  `variable_config`.

- anchor_draw:

  How anchors fire. `"all"` (default) fires every anchor, minus
  `anchor_dropout`. `"prevalence"` fires each anchor independently with
  probability `prevalence * (1 - anchor_dropout)`, so anchors at
  prevalence 0.8 are absent from about 20% of images, as in real
  samples. Needs a `prevalence` column.

## Value

List with `clean` (H x W), `noisy` (H x W), and per-layer composition
metadata (`anchor_ids`: the anchors that fired this sample, after any
`anchor_draw` / `anchor_dropout` thinning; `variable_ids`,
`contamination_positions`, and `rendered_ids` — the compound_ids
actually drawn in the anchors / variable layers, used as ground-truth
presence labels).

## Details

1.  **Anchor backbone.** Every compound in `anchor_ids` fires (subject
    to `anchor_dropout` / `anchor_draw`), with its own catalog
    morphology and normal per-peak location / size jitter. These are the
    stable ~100 highly-prevalent compounds you want the encoder to treat
    as the baseline landscape.

2.  **Variable catalog compounds** (optional). Per-sample random subset
    of the catalog's less-prevalent compounds — filtered by
    `variable_config$prevalence_range` and (optionally) excluding the
    anchor pool. Simulates diet-related, host-genetics-related, and
    other per-person volatiles. Intensity scaled by a per-peak factor
    drawn from `variable_config$intensity_scale`. The subset is drawn
    uniformly (every eligible compound equally likely) or, with
    `variable_config$draw = "prevalence"`, each compound fires at its
    own catalog prevalence.

3.  **Contamination** (optional). Uniform-random peaks placed at
    positions NOT overlapping catalog compounds. Simulates sampling
    apparatus contamination, environmental exposure, and truly
    unfamiliar peaks. Sigmas drawn from the pooled catalog sigma
    distribution so shapes look realistic.
