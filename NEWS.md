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

## Performance

* `pretrain_denoising_online(num_workers > 1)` no longer re-sends its state
  to the workers every batch. Each batch used to export closures defined
  inside the function, which serialise the whole function frame (catalog,
  `val_real`, model, optimiser; under `do.call()` the catalog and
  `val_real` twice). That came to ~0.2 GiB per future with 2 `val_real`
  samples and ~0.5 GiB with 30, above future's default 500 MiB
  `future.globals.maxSize`. Workers are now persistent PSOCK processes
  that receive the catalog and layer settings once per run; each batch
  sends only the per-sample augmentation specs and seeds (KB to a few
  MB). On the 1734-compound universal catalog (1504 x 96, batch 32,
  4 workers, MPS) peak resident memory of main + workers fell from
  7.2 GB to 4.2 GB and worker memory no longer saw-tooths.
* Peak rendering in `generate_one_synthetic_from_catalog()` and
  `generate_one_synthetic()` adds each peak in place instead of copying
  the whole H x W canvas per peak: ~4x faster per image (0.43 s to
  0.10 s on the realistic recipe) with identical output.
* Together: 6 steps of the realistic recipe at batch 32 with 4 workers
  went from 75 s to 19 s.

## Behaviour changes

* `pretrain_denoising_online(num_workers > 1)` no longer sets (or
  restores) the global `future::plan()`; it dispatches with
  `parallel::clusterApply()`. The per-sample L'Ecuyer-CMRG streams are
  the ones `future_lapply(future.seed = TRUE)` assigned, so batches are
  bit-for-bit the same as before for a fixed `seed`, serial and
  parallel. `future.apply` moved from Imports to Suggests.
* Parallel workers run the installed classydms; a warning now says so
  when it differs from the code loaded in the session (e.g. under
  `devtools::load_all()`).

## Documentation

* `contamination_config` docs now note that the default
  `intensity_range = c(0.1, 0.5)` and exclusion boxes are much brighter
  and wider than real off-catalog peaks (median normalised intensity
  about 0.014 on a 606-sample universal catalog).

## Bug fixes

* `pretrain_denoising_online()` with `contamination_schedule` no longer
  fails ("non-finite pixels") when `contamination_config` omits
  `n_per_sample` or `intensity_range`. The light regime scaled the
  missing value (`NULL * 0.3` is `numeric(0)`), so the generator's
  defaults never applied; they are now filled in before scaling.
* `generate_noisy_pretrain_sample()`: `variable_config$exclude_anchors`
  now excludes every compound in `anchor_ids`, not just the anchors that
  survived `anchor_draw = "prevalence"` / `anchor_dropout` this sample.
  Before, a dropped anchor whose prevalence fell inside
  `variable_config$prevalence_range` got a second chance in the variable
  layer, so with both draws on `"prevalence"` an anchor of prevalence `p`
  fired with probability `1 - (1 - p)^2` instead of `p` (uniform mode
  with `anchor_dropout > 0` leaked too, by less). **Output changes** for
  the same seed whenever anchors are thinned and the variable range
  covers them; output with no thinning, or with a variable range below
  every anchor, is bit-for-bit unchanged. The returned `anchor_ids` are
  still the anchors that fired.

* The uniform variable draw no longer misbehaves when exactly one
  compound is eligible (`sample()` on a length-1 vector sampled from
  `1:id` instead of returning the id).
