# classydms (development version)

## New features

* `generate_noisy_pretrain_sample()` and `pretrain_denoising_online()`
  gain a prevalence-weighted draw for the variable layer:
  `variable_config$draw = "prevalence"` fires each eligible compound
  independently with probability `min(1, prevalence * prevalence_scale)`
  (`prevalence_scale` defaults to 1), so per-compound frequencies and the
  expected count per image follow the catalog. In this mode
  `n_per_sample` is an optional per-sample cap (default: no cap). The
  default `draw = "uniform"` reproduces earlier output exactly for the
  same seed.
* New `anchor_draw = "prevalence"` argument fires each anchor with
  probability `prevalence * (1 - anchor_dropout)` instead of always
  (default `anchor_draw = "all"`, unchanged behaviour).
* `pretrain_denoising_online()` checks the draw settings at entry and
  records `anchor_draw`, `variable_draw` and `variable_prevalence_scale`
  in the manifest `hyperparams`.

## Documentation

* `contamination_config` docs now note that the default
  `intensity_range = c(0.1, 0.5)` and exclusion boxes are much brighter
  and wider than real off-catalog peaks (median normalised intensity
  about 0.014 on a 606-sample universal catalog).

## Bug fixes

* The uniform variable draw no longer misbehaves when exactly one
  compound is eligible (`sample()` on a length-1 vector sampled from
  `1:id` instead of returning the id).
