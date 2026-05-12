# Reconstruct real images through a pre-trained autoencoder

Feeds each requested sample through the autoencoder and returns the
original, the reconstruction, and the per-sample MSE. Useful as a
qualitative check that synthetic pretraining produced an autoencoder
that generalizes to real data.

## Usage

``` r
reconstruct_real_samples(
  autoencoder,
  Z_list,
  sample_indices = 1:6,
  device = "cpu"
)
```

## Arguments

- autoencoder:

  A trained `dms_denoising_autoencoder`.

- Z_list:

  List of preprocessed, padded intensity matrices.

- sample_indices:

  Which samples to reconstruct (default first 6).

- device:

  "cpu" or "cuda".

## Value

List of `list(index, original, reconstructed, mse)`.
