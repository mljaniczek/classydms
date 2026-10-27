# Pre-train the encoder from a peak catalog (smoke test only)

Same idea as
[`pretrain_autoencoder()`](https://mljaniczek.github.io/classydms/reference/pretrain_autoencoder.md)
but generates the synthetic pretraining set via
[`generate_one_synthetic_from_catalog()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic_from_catalog.md).
Use this when you have a catalog (perhaps loaded from disk) but not the
original `peak_params` object.

## Usage

``` r
pretrain_autoencoder_from_catalog(
  catalog,
  H,
  W,
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

- n_synthetic:

  Number of synthetic samples to generate.

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

List with `encoder`, `autoencoder`, `loss_history`.

## Details

**For real pretraining runs prefer
[`pretrain_denoising_online()`](https://mljaniczek.github.io/classydms/reference/pretrain_denoising_online.md)
with `catalog = your_catalog`** — that generates fresh batches per
training step, so the model can see effectively unlimited unique
synthetic samples with constant memory. This pre-alloc version holds the
full `n_synthetic` × H × W tensor in RAM and reuses the same set every
epoch, which caps the effective training data. It stays available for
smoke tests (`n_synthetic <= a few thousand`).
