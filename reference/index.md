# Package index

## Loading data

Read raw GC-DMS files and assemble a samples list.

- [`read_dms_file()`](https://mljaniczek.github.io/classydms/reference/read_dms_file.md)
  : Read one raw GC-DMS data file
- [`load_dms_directory()`](https://mljaniczek.github.io/classydms/reference/load_dms_directory.md)
  : Read all GC-DMS files in a directory and assemble a samples list

## Preprocessing

Smoothing, baseline correction, dust thresholding, trimming,
normalization, and padding for GC-DMS intensity matrices.

- [`preprocess_matrix_rt()`](https://mljaniczek.github.io/classydms/reference/preprocess_matrix_rt.md)
  : Savitzky-Golay smoothing + AsLS baseline correction (per CV column)

- [`process_one_sample()`](https://mljaniczek.github.io/classydms/reference/process_one_sample.md)
  :

  Run preprocessing on one sample (list with `$Z`, optional `$time`,
  `$cv`)

- [`baseline_basement()`](https://mljaniczek.github.io/classydms/reference/baseline_basement.md)
  : Dust-threshold a sample below a fixed intensity

- [`trim_bounds_from_occupancy()`](https://mljaniczek.github.io/classydms/reference/trim_bounds_from_occupancy.md)
  : Compute occupancy-based trim bounds for one sample

- [`trim_sample()`](https://mljaniczek.github.io/classydms/reference/trim_sample.md)
  : Apply trim bounds to a sample (clamps to actual size if smaller)

- [`pooled_trim_bounds()`](https://mljaniczek.github.io/classydms/reference/pooled_trim_bounds.md)
  : Pool per-sample trim bounds to a single cohort-wide bounding box

- [`normalize_sample()`](https://mljaniczek.github.io/classydms/reference/normalize_sample.md)
  : Log-quantile normalize a single sample

- [`compute_pad_targets()`](https://mljaniczek.github.io/classydms/reference/compute_pad_targets.md)
  : Compute padding targets across a list of matrices

- [`pad_to_target()`](https://mljaniczek.github.io/classydms/reference/pad_to_target.md)
  : Center-pad a matrix to target dimensions

- [`robust_eps()`](https://mljaniczek.github.io/classydms/reference/robust_eps.md)
  : Compute a per-sample relative "dust" threshold

## Synthetic data generation

Estimate peak parameter distributions from real data (label-agnostic)
and generate synthetic GC-DMS images statistically matched to the real
cohort.

- [`estimate_peak_params()`](https://mljaniczek.github.io/classydms/reference/estimate_peak_params.md)
  : Estimate peak parameter distributions from real GC-DMS samples
- [`generate_one_synthetic()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic.md)
  : Generate one synthetic GC-DMS image
- [`generate_synthetic_dataset()`](https://mljaniczek.github.io/classydms/reference/generate_synthetic_dataset.md)
  : Generate N synthetic pairs as torch tensors

## Autoencoder architecture

ResNet-style encoder, mirrored decoder, and composed autoencoder.

- [`dms_encoder()`](https://mljaniczek.github.io/classydms/reference/dms_encoder.md)
  : GC-DMS encoder (ResNet-style, asymmetric stride)
- [`dms_decoder()`](https://mljaniczek.github.io/classydms/reference/dms_decoder.md)
  : GC-DMS decoder (mirror of encoder, transposed convolutions)
- [`dms_denoising_autoencoder()`](https://mljaniczek.github.io/classydms/reference/dms_denoising_autoencoder.md)
  : Composed denoising autoencoder (encoder + decoder)
- [`res_block()`](https://mljaniczek.github.io/classydms/reference/res_block.md)
  : Residual block used by encoder and tiny ResNet variants

## Pre-training

Self-supervised training of the encoder on synthetic data, either via
the in-memory batch trainer or the online variant for scale-up.

- [`pretrain_autoencoder()`](https://mljaniczek.github.io/classydms/reference/pretrain_autoencoder.md)
  [`pretrain_denoising()`](https://mljaniczek.github.io/classydms/reference/pretrain_autoencoder.md)
  : Pre-train the encoder via batch denoising autoencoder
- [`pretrain_denoising_online()`](https://mljaniczek.github.io/classydms/reference/pretrain_denoising_online.md)
  : Pre-train the encoder via online denoising autoencoder

## Classification

Use the pre-trained encoder as a frozen feature extractor for downstream
classification on real labeled samples.

- [`build_pretrained_classifier()`](https://mljaniczek.github.io/classydms/reference/build_pretrained_classifier.md)
  : Attach a classifier head to a pre-trained encoder
- [`train_one_fold_pretrained()`](https://mljaniczek.github.io/classydms/reference/train_one_fold_pretrained.md)
  : Train one classifier fold with a pre-trained encoder
- [`cv_pretrained_oof()`](https://mljaniczek.github.io/classydms/reference/cv_pretrained_oof.md)
  : 5-fold CV classification with a pre-trained encoder
- [`extract_encoder_features()`](https://mljaniczek.github.io/classydms/reference/extract_encoder_features.md)
  : Extract 64-dim encoder features for all samples (GAP-pooled)
- [`cv_encoder_elastic_net()`](https://mljaniczek.github.io/classydms/reference/cv_encoder_elastic_net.md)
  : Elastic net classification on 64-dim encoder features

## Interpretability and QC

Reconstruction QC, spatial saliency mapping, and aggregation tools for
tracing classifier decisions back to (RT, CV) coordinates.

- [`reconstruct_real_samples()`](https://mljaniczek.github.io/classydms/reference/reconstruct_real_samples.md)
  : Reconstruct real images through a pre-trained autoencoder
- [`compute_saliency_map()`](https://mljaniczek.github.io/classydms/reference/compute_saliency_map.md)
  : Compute a Class Activation Map (CAM) for one sample
- [`aggregate_saliency_masked()`](https://mljaniczek.github.io/classydms/reference/aggregate_saliency_masked.md)
  : Aggregate saliency maps across samples, masking padded regions

## Pre-trained model helpers

Download and cache pre-trained encoders distributed via GitHub Releases.

- [`load_pretrained_encoder()`](https://mljaniczek.github.io/classydms/reference/load_pretrained_encoder.md)
  : Load a pre-trained encoder by name
- [`list_pretrained_encoders()`](https://mljaniczek.github.io/classydms/reference/list_pretrained_encoders.md)
  : List available pre-trained encoder files
