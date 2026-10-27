# Changelog

## classydms (development version)

### New features

- [`generate_noisy_pretrain_sample()`](https://mljaniczek.github.io/classydms/reference/generate_noisy_pretrain_sample.md)
  and
  [`pretrain_denoising_online()`](https://mljaniczek.github.io/classydms/reference/pretrain_denoising_online.md)
  gain a prevalence-weighted draw for the variable layer:
  `variable_config$draw = "prevalence"` fires each eligible compound
  independently with probability `min(1, prevalence * prevalence_scale)`
  (`prevalence_scale` defaults to 1), so per-compound frequencies and
  the expected count per image follow the catalog. In this mode
  `n_per_sample` is an optional per-sample cap (default: no cap). The
  default `draw = "uniform"` reproduces earlier output exactly for the
  same seed.
- New `anchor_draw = "prevalence"` argument fires each anchor with
  probability `prevalence * (1 - anchor_dropout)` instead of always
  (default `anchor_draw = "all"`, unchanged behaviour).
- [`pretrain_denoising_online()`](https://mljaniczek.github.io/classydms/reference/pretrain_denoising_online.md)
  checks the draw settings at entry and records `anchor_draw`,
  `variable_draw` and `variable_prevalence_scale` in the manifest
  `hyperparams`.

### Documentation

- `contamination_config` docs now note that the default
  `intensity_range = c(0.1, 0.5)` and exclusion boxes are much brighter
  and wider than real off-catalog peaks (median normalised intensity
  about 0.014 on a 606-sample universal catalog).

### Bug fixes

- [`generate_noisy_pretrain_sample()`](https://mljaniczek.github.io/classydms/reference/generate_noisy_pretrain_sample.md):
  `variable_config$exclude_anchors` now excludes every compound in
  `anchor_ids`, not just the anchors that survived
  `anchor_draw = "prevalence"` / `anchor_dropout` this sample. Before, a
  dropped anchor whose prevalence fell inside
  `variable_config$prevalence_range` got a second chance in the variable
  layer, so with both draws on `"prevalence"` an anchor of prevalence
  `p` fired with probability `1 - (1 - p)^2` instead of `p` (uniform
  mode with `anchor_dropout > 0` leaked too, by less). **Output
  changes** for the same seed whenever anchors are thinned and the
  variable range covers them; output with no thinning, or with a
  variable range below every anchor, is bit-for-bit unchanged. The
  returned `anchor_ids` are still the anchors that fired.

- The uniform variable draw no longer misbehaves when exactly one
  compound is eligible ([`sample()`](https://rdrr.io/r/base/sample.html)
  on a length-1 vector sampled from `1:id` instead of returning the id).
