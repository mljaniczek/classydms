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

## Quality control

Per-sample QC features and exclusion criteria, applied between
preprocessing and peak parameter estimation.

- [`qc_features()`](https://mljaniczek.github.io/classydms/reference/qc_features.md)
  : Compute per-sample QC features
- [`flag_junk_samples()`](https://mljaniczek.github.io/classydms/reference/flag_junk_samples.md)
  : Flag samples that fail QC checks
- [`plot_qc()`](https://mljaniczek.github.io/classydms/reference/plot_qc.md)
  : Visualize per-sample QC features

## Common physical grid

Put every cohort on one (seconds, volts) grid so position never encodes
cohort.

- [`dms_grid()`](https://mljaniczek.github.io/classydms/reference/dms_grid.md)
  : Define a common physical grid for GC-DMS images
- [`make_dms_grid()`](https://mljaniczek.github.io/classydms/reference/make_dms_grid.md)
  : Build a dms_grid from start, step and size
- [`choose_cv_window()`](https://mljaniczek.github.io/classydms/reference/choose_cv_window.md)
  : Choose a common CV window from pooled signal occupancy
- [`resample_to_grid()`](https://mljaniczek.github.io/classydms/reference/resample_to_grid.md)
  : Resample one sample onto a dms_grid (RT and CV)
- [`resample_to_cv_grid()`](https://mljaniczek.github.io/classydms/reference/resample_to_cv_grid.md)
  : Resample one sample onto a CV grid
- [`resample_to_time_grid()`](https://mljaniczek.github.io/classydms/reference/resample_to_time_grid.md)
  : Resample one sample onto a regular time grid via linear
  interpolation
- [`compute_common_time_grid()`](https://mljaniczek.github.io/classydms/reference/compute_common_time_grid.md)
  : Compute a common regular time grid across a cohort of samples

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
- [`synthetic_quality_check()`](https://mljaniczek.github.io/classydms/reference/synthetic_quality_check.md)
  : Quantitative comparison of synthetic vs real GC-DMS samples
- [`plot_synthetic_quality()`](https://mljaniczek.github.io/classydms/reference/plot_synthetic_quality.md)
  : Plot real vs synthetic per-sample summary statistics

## Peak parameters and alignment

Validate, subset, merge and align peak_params bundles across samples and
cohorts.

- [`validate_peak_params()`](https://mljaniczek.github.io/classydms/reference/validate_peak_params.md)
  : Validate the internal consistency of a peak_params object
- [`subset_peak_params()`](https://mljaniczek.github.io/classydms/reference/subset_peak_params.md)
  : Subset a peak_params object to a chosen set of samples
- [`merge_peak_params()`](https://mljaniczek.github.io/classydms/reference/merge_peak_params.md)
  : Merge multiple peak_params objects into a single super-cohort
- [`backfill_sample_idx_raw()`](https://mljaniczek.github.io/classydms/reference/backfill_sample_idx_raw.md)
  : Backfill sample_idx_raw and sample_names on an older peak_params
- [`align_across_cohorts()`](https://mljaniczek.github.io/classydms/reference/align_across_cohorts.md)
  : Align samples across multiple cohorts jointly
- [`align_via_catalog()`](https://mljaniczek.github.io/classydms/reference/align_via_catalog.md)
  : Align samples via cluster-membership offsets
- [`render_peak_params_report()`](https://mljaniczek.github.io/classydms/reference/render_peak_params_report.md)
  : Render the peak-parameters diagnostic report

## Peak catalogs

Cluster peaks into shared compounds, inspect catalogs, and generate
synthetic samples from them.

- [`build_peak_catalog()`](https://mljaniczek.github.io/classydms/reference/build_peak_catalog.md)
  : Build a peak catalog from a peak_params object
- [`build_universal_catalog()`](https://mljaniczek.github.io/classydms/reference/build_universal_catalog.md)
  : Build a universal catalog by pooling all samples from multiple
  cohorts
- [`merge_catalogs()`](https://mljaniczek.github.io/classydms/reference/merge_catalogs.md)
  : Merge multiple peak catalogs into a single unified catalog
- [`split_catalog()`](https://mljaniczek.github.io/classydms/reference/split_catalog.md)
  : Split a catalog into two disjoint sub-catalogs by compound
- [`compare_catalogs()`](https://mljaniczek.github.io/classydms/reference/compare_catalogs.md)
  : Compare two peak catalogs by compound proximity
- [`catalog_summary()`](https://mljaniczek.github.io/classydms/reference/catalog_summary.md)
  : Diagnostic summary of a peak catalog
- [`catalog_compactness()`](https://mljaniczek.github.io/classydms/reference/catalog_compactness.md)
  : Compactness diagnostic for a peak catalog
- [`catalog_pixel_positions()`](https://mljaniczek.github.io/classydms/reference/catalog_pixel_positions.md)
  : Catalog compound positions on a padded canvas
- [`plot_catalog()`](https://mljaniczek.github.io/classydms/reference/plot_catalog.md)
  : Plot the spatial distribution of a peak catalog
- [`generate_one_synthetic_from_catalog()`](https://mljaniczek.github.io/classydms/reference/generate_one_synthetic_from_catalog.md)
  : Generate one synthetic sample from a peak catalog
- [`generate_noisy_pretrain_sample()`](https://mljaniczek.github.io/classydms/reference/generate_noisy_pretrain_sample.md)
  : Generate one three-layer synthetic sample for encoder pretraining

## Case/control simulation

Simulate labelled cohorts with known injected biomarkers to test whether
the pipeline recovers them.

- [`simulate_case_control_cohort()`](https://mljaniczek.github.io/classydms/reference/simulate_case_control_cohort.md)
  : Simulate a case / control cohort with known injected biomarkers
- [`simulate_case_control_from_catalog()`](https://mljaniczek.github.io/classydms/reference/simulate_case_control_from_catalog.md)
  : Simulate a case / control cohort from a peak catalog
- [`holdout_biomarkers()`](https://mljaniczek.github.io/classydms/reference/holdout_biomarkers.md)
  : Hold out K observed peaks as biomarker templates
- [`build_biomarker_spec()`](https://mljaniczek.github.io/classydms/reference/build_biomarker_spec.md)
  : Combine biomarker templates with case/control differential spec
- [`build_biomarker_spec_catalog()`](https://mljaniczek.github.io/classydms/reference/build_biomarker_spec_catalog.md)
  : Build a case/control biomarker spec for catalog-based simulation
- [`pick_catalog_biomarkers()`](https://mljaniczek.github.io/classydms/reference/pick_catalog_biomarkers.md)
  : Pick a set of catalog compounds as biomarkers
- [`biomarkers_level1()`](https://mljaniczek.github.io/classydms/reference/biomarkers_level1.md)
  : Level 1 (sanity check) biomarker scenario
- [`biomarkers_level2()`](https://mljaniczek.github.io/classydms/reference/biomarkers_level2.md)
  : Level 2 (realistic) biomarker scenario
- [`biomarkers_level3()`](https://mljaniczek.github.io/classydms/reference/biomarkers_level3.md)
  : Level 3 (hard) biomarker scenario
- [`biomarkers_sanity_max()`](https://mljaniczek.github.io/classydms/reference/biomarkers_sanity_max.md)
  : Maximally-favorable ("sanity-max") biomarker scenario
- [`catalog_biomarkers_sanity_max()`](https://mljaniczek.github.io/classydms/reference/catalog_biomarkers_sanity_max.md)
  : Sanity-max scenario for catalog-based case/control simulation
- [`catalog_biomarkers_universal_backdrop()`](https://mljaniczek.github.io/classydms/reference/catalog_biomarkers_universal_backdrop.md)
  : Universal-backdrop biomarker preset for catalog simulation
- [`evaluate_biomarker_recovery()`](https://mljaniczek.github.io/classydms/reference/evaluate_biomarker_recovery.md)
  : Evaluate biomarker recovery from a saliency map

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
- [`pretrain_autoencoder_from_catalog()`](https://mljaniczek.github.io/classydms/reference/pretrain_autoencoder_from_catalog.md)
  : Pre-train the encoder from a peak catalog (smoke test only)
- [`pretrain_autoencoder_from_catalog_noisy()`](https://mljaniczek.github.io/classydms/reference/pretrain_autoencoder_from_catalog_noisy.md)
  : Pre-train the encoder from a peak catalog using three-layer noisy
  samples (smoke test only)

## Compound-presence head and readouts

Catalog-anchored attention pooling and the presence head trained
alongside the autoencoder.

- [`dms_presence_head()`](https://mljaniczek.github.io/classydms/reference/dms_presence_head.md)
  : Compound-presence head
- [`compound_head_geometry()`](https://mljaniczek.github.io/classydms/reference/compound_head_geometry.md)
  : Catalog geometry for a compound-anchored head
- [`compound_attention_pool()`](https://mljaniczek.github.io/classydms/reference/compound_attention_pool.md)
  : Compound-anchored attention pooling
- [`attention_pool_readout()`](https://mljaniczek.github.io/classydms/reference/attention_pool_readout.md)
  : Learned-query attention pooling readout

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
  : Extract encoder features for all samples
- [`cv_encoder_elastic_net()`](https://mljaniczek.github.io/classydms/reference/cv_encoder_elastic_net.md)
  : Elastic net classification on encoder features

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
