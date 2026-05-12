# classydms

GC-DMS Disease Classification via Synthetic Pre-Training

`classydms` is an R package that implements an end-to-end pipeline for disease vs control classification from Gas Chromatography Differential Mobility Spectrometry (GC-DMS) data. It is designed for small-sample clinical studies (n ~ 100-200 subjects) and supports applying the same workflow across many disease settings.

## Pipeline overview

```
Raw GC-DMS file (RT x CV intensity matrix)
        |
        v
Preprocessing: Savitzky-Golay + AsLS + dust threshold + occupancy trim + log-quantile normalize + zero-pad
        |
        v
Estimate peak parameter distributions (label-agnostic)
        |
        v
Generate millions of synthetic GC-DMS images matching the distribution
        |
        v
Pre-train a ResNet-style encoder via self-supervised denoising autoencoder
        |
        v
Extract 64-dim feature vector per real sample (frozen encoder + global average pool)
        |
        v
Elastic net classifier on the 64 features  -->  cross-validated AUC, interpretable coefficients
        |
        v
Spatial saliency mapping (CAM weighted by elastic net coefficients) --> RT x CV biomarker hotspots
```

## Installation

```r
# Install development version from GitHub
# install.packages("remotes")
remotes::install_github("mljaniczek/classydms")
```

The package depends on `torch` (R bindings to LibTorch). If this is your first time installing `torch`, run `torch::install_torch()` after package installation.

## Quick start

```r
library(classydms)

# 1. Load raw data (samples_keep, y_keep) using your project-specific loader
# 2. Preprocess
samples_proc <- lapply(samples_keep, process_one_sample)
b_list <- lapply(samples_proc, function(s) trim_bounds_from_occupancy(s$Z))
R0 <- max(sapply(samples_proc, function(s) nrow(s$Z)))
C0 <- max(sapply(samples_proc, function(s) ncol(s$Z)))
global_b <- pooled_trim_bounds(b_list, R = R0, C = C0)
samples_trimmed <- lapply(samples_proc, trim_sample, bounds = global_b)

Z_norm   <- purrr::map(samples_trimmed, ~normalize_sample(.x$Z))
pad_dims <- compute_pad_targets(Z_norm)
Z_padded <- purrr::map(Z_norm, ~pad_to_target(.x, pad_dims$H, pad_dims$W))

# 3. Estimate peak parameters from real data (label-agnostic)
Z_raw     <- purrr::map(samples_trimmed, ~.x$Z)
Z_pretrim <- purrr::map(samples_trimmed,
                          ~if (!is.null(.x$Z_pretrim)) .x$Z_pretrim else .x$Z)
peak_params <- estimate_peak_params(Z_raw, Z_pretrim_list = Z_pretrim,
                                     top_k = 150, frac_height = 0.5)

# 4. Self-supervised pretraining (recommended: online, for memory efficiency)
pretrain_result <- pretrain_denoising_online(
  peak_params,
  H = pad_dims$H, W = pad_dims$W,
  steps_per_epoch = 1000L, epochs = 50L, batch_size = 32L,
  save_path = "encoder.pt"
)

# 5. Classify with elastic net on encoder features
elnet <- cv_encoder_elastic_net(pretrain_result$encoder,
                                 Z_list = Z_padded, y = y_keep,
                                 k = 5, alpha = 0.5)
elnet$auc
elnet$coefs    # which of the 64 channels matter

# 6. Spatial saliency map for biomarker localization
disease_idx <- which(y_keep == 2)
orig_dims <- lapply(samples_trimmed, function(s) c(H = nrow(s$Z), W = ncol(s$Z)))
sal <- aggregate_saliency_masked(pretrain_result$encoder, Z_padded,
                                   coefs = elnet$coefs,
                                   sample_indices = disease_idx,
                                   orig_dims = orig_dims)
# sal is a 2-D matrix in padded image coordinates; map back to RT/CV using
# samples_trimmed[[1]]$time and samples_trimmed[[1]]$cv
```

## Using a pre-trained encoder

To save weeks of compute, you can download an encoder trained on a previous disease cohort:

```r
# List available encoders attached to GitHub Releases
list_pretrained_encoders()

# Download and use a specific encoder
encoder <- load_pretrained_encoder("flu_v1")

# Or fetch the full autoencoder for reconstruction QC
bundle <- load_pretrained_encoder("flu_v1", with_autoencoder = TRUE)
bundle$manifest$hyperparams      # how it was trained
bundle$manifest$peak_params      # what synthetic distribution it saw
```

Encoders are cached locally on first download.

## Documentation

- Vignettes: see `vignette(package = "classydms")` for end-to-end walk-throughs.
- Function reference: `?pretrain_denoising_online`, `?cv_encoder_elastic_net`, `?compute_saliency_map`, etc.
- Methods: full preprocessing rationale and pipeline description in the proposal document.

## Citation

If you use `classydms` in your work, please cite (TODO: add citation once published).

## License

MIT
