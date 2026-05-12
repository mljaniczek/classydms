# Pre-train the encoder via batch denoising autoencoder

Generates `n_synthetic` synthetic GC-DMS images, trains a denoising
autoencoder on them (noisy -\> clean), and returns the trained encoder.
Memory-limited: the full synthetic dataset is held in RAM. For larger
totals, use
[`pretrain_denoising_online()`](https://mljaniczek.github.io/classydms/reference/pretrain_denoising_online.md).

## Usage

``` r
pretrain_autoencoder(
  peak_params,
  H,
  W,
  n_synthetic = 5000L,
  epochs = 50L,
  batch_size = 32L,
  lr = 0.001,
  weight_decay = 1e-04,
  add_noise = TRUE,
  stem_stride_rt = 4L,
  stem_stride_cv = 1L,
  device = if (torch::cuda_is_available()) "cuda" else "cpu",
  save_path = NULL,
  seed = 42L
)

pretrain_denoising(
  peak_params,
  H,
  W,
  n_synthetic = 5000L,
  epochs = 50L,
  batch_size = 32L,
  lr = 0.001,
  weight_decay = 1e-04,
  add_noise = TRUE,
  stem_stride_rt = 4L,
  stem_stride_cv = 1L,
  device = if (torch::cuda_is_available()) "cuda" else "cpu",
  save_path = NULL,
  seed = 42L
)
```

## Arguments

- peak_params:

  Output of
  [`estimate_peak_params()`](https://mljaniczek.github.io/classydms/reference/estimate_peak_params.md).

- H, W:

  Synthetic image dimensions (must match downstream input).

- n_synthetic:

  Number of synthetic images.

- epochs, batch_size, lr, weight_decay:

  Optimization parameters.

- add_noise:

  If TRUE, denoising AE; if FALSE, reconstruction AE.

- stem_stride_rt, stem_stride_cv:

  Stem stride for the encoder.

- device:

  "cpu" or "cuda".

- save_path:

  If non-NULL, writes encoder and full autoencoder to disk.

- seed:

  RNG seed.

## Value

List with `encoder`, `autoencoder`, and `loss_history`.

## See also

[`pretrain_denoising_online()`](https://mljaniczek.github.io/classydms/reference/pretrain_denoising_online.md)
