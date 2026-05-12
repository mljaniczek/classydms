# Pre-train the encoder via online denoising autoencoder

Generates fresh synthetic batches on-the-fly during training, giving
constant memory usage regardless of total samples seen. Recommended for
scaling to hundreds of thousands or millions of total samples.

## Usage

``` r
pretrain_denoising_online(
  peak_params,
  H,
  W,
  steps_per_epoch = 500L,
  epochs = 30L,
  batch_size = 32L,
  lr = 0.001,
  weight_decay = 1e-04,
  add_noise = TRUE,
  size_jitter = 0.6,
  grad_clip = 1,
  norm_clamp = 10,
  val_n = 200L,
  checkpoint_every = 10L,
  loss_diagnostics = TRUE,
  stem_stride_rt = 4L,
  stem_stride_cv = 1L,
  device = if (torch::cuda_is_available()) "cuda" else "cpu",
  save_path = NULL,
  resume_from = NULL,
  seed = 42L,
  verbose = TRUE
)
```

## Arguments

- peak_params:

  Output of
  [`estimate_peak_params()`](https://mljaniczek.github.io/classydms/reference/estimate_peak_params.md).

- H, W:

  Synthetic image dimensions.

- steps_per_epoch:

  Number of training batches per epoch.

- epochs:

  Number of epochs.

- batch_size:

  Batch size.

- lr, weight_decay:

  Optimizer parameters.

- add_noise:

  If TRUE, denoising AE; if FALSE, reconstruction AE.

- size_jitter:

  Per-peak random log-normal size scale factor SD.

- grad_clip:

  Maximum L2 norm for gradient clipping (0 disables).

- norm_clamp:

  Maximum normalized pixel value (safety clamp).

- val_n:

  Number of fixed validation samples (0 disables validation).

- checkpoint_every:

  Save encoder checkpoint every N epochs.

- loss_diagnostics:

  If TRUE, also print median, p95, max batch loss.

- stem_stride_rt, stem_stride_cv:

  Encoder stem strides.

- device:

  "cpu" or "cuda".

- save_path:

  If non-NULL, writes encoder, autoencoder, and manifest to disk after
  every epoch (so a mid-run crash leaves recoverable artifacts).

- resume_from:

  If non-NULL, path to a previously-saved encoder (`*.pt`). The function
  looks for the corresponding `_autoencoder.pt` and `_manifest.Rdata`
  next to it and resumes training from `last_epoch_completed + 1`,
  inheriting prior loss history. Useful for picking up after a crash
  without losing the work done so far. The other hyperparameters (`H`,
  `W`, `epochs`, etc.) must match the original run.

- seed:

  RNG seed.

- verbose:

  Whether to print epoch progress.

## Value

List with `encoder`, `autoencoder`, `loss_history`, `val_loss_history`,
`batch_loss_stats`, and `total_samples`.

## Details

Includes scale-up safety features: a fixed validation set held constant
across epochs for stable progress tracking, periodic encoder
checkpointing for crash recovery, per-epoch manifest save with all
hyperparameters for reproducibility, and per-epoch ETA estimate.

Total samples seen = `steps_per_epoch * epochs * batch_size`.
