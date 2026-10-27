# Pre-train the encoder from a peak catalog using three-layer noisy samples (smoke test only)

Same training loop as
[`pretrain_autoencoder_from_catalog()`](https://mljaniczek.github.io/classydms/reference/pretrain_autoencoder_from_catalog.md)
but each synthetic sample is generated via
[`generate_noisy_pretrain_sample()`](https://mljaniczek.github.io/classydms/reference/generate_noisy_pretrain_sample.md)
— an anchor backbone plus optional variable-catalog and contamination
layers. The extra variability is intended to teach the encoder to expect
real-world skin-VOC characteristics: consistent backbone, per-sample
dietary / host-genetic peaks, and contamination-like unfamiliar peaks.

## Usage

``` r
pretrain_autoencoder_from_catalog_noisy(
  catalog,
  H,
  W,
  anchor_ids,
  variable_config = NULL,
  contamination_config = NULL,
  n_synthetic = 2000L,
  epochs = 30L,
  batch_size = 32L,
  lr = 0.001,
  weight_decay = 1e-04,
  add_noise = TRUE,
  stem_stride_rt = 4L,
  stem_stride_cv = 1L,
  device = if (torch::cuda_is_available()) "cuda" else "cpu",
  seed = 42L
)
```

## Arguments

- catalog:

  A `peak_catalog` object.

- H, W:

  Synthetic image dimensions.

- anchor_ids:

  Integer vector of compound_ids to force-fire in every sample (the
  invariant backbone). Typically the ~100 anchors at prevalence \>= 0.9
  from the universal catalog.

- variable_config, contamination_config:

  Passed through to
  [`generate_noisy_pretrain_sample()`](https://mljaniczek.github.io/classydms/reference/generate_noisy_pretrain_sample.md).
  Pass `NULL` for either to skip that layer. See that function's docs
  for field details.

- n_synthetic:

  Number of synthetic samples generated.

- epochs, batch_size, lr, weight_decay:

  Optimization parameters.

- add_noise:

  If TRUE, denoising autoencoder; if FALSE, reconstruction.

- stem_stride_rt, stem_stride_cv:

  Encoder stem strides.

- device:

  "cpu" or "cuda".

- seed:

  RNG seed.

## Value

List with `encoder`, `autoencoder`, `loss_history`, and the
`variable_config` / `contamination_config` actually used (echoed for
provenance).

## Details

**For real pretraining runs prefer
[`pretrain_denoising_online()`](https://mljaniczek.github.io/classydms/reference/pretrain_denoising_online.md)
with `catalog = your_catalog`, `anchor_ids = ids`, and the same
`variable_config` / `contamination_config`** — that path is online
(fresh batches per step) so training data isn't capped at `n_synthetic`.
This pre-alloc version stays available for smoke tests.

To reproduce the current
[`pretrain_autoencoder_from_catalog()`](https://mljaniczek.github.io/classydms/reference/pretrain_autoencoder_from_catalog.md)
behavior, pass `variable_config = NULL` and
`contamination_config = NULL` — the samples will contain only the anchor
set (fired at 100 percent).
