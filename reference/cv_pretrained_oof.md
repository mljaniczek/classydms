# 5-fold CV classification with a pre-trained encoder

Out-of-fold predictions are aggregated into a single AUC. Each fold
independently deep-clones the encoder via
[`build_pretrained_classifier()`](https://mljaniczek.github.io/classydms/reference/build_pretrained_classifier.md)
to prevent cross-fold leakage.

## Usage

``` r
cv_pretrained_oof(
  Z_list,
  y,
  encoder,
  k = 5L,
  epochs = 20L,
  batch_size = 4L,
  lr = 0.001,
  weight_decay = 1e-04,
  dropout_p = 0.3,
  freeze_backbone = TRUE,
  unfreeze_last_block = TRUE,
  augment = TRUE,
  rt_shift = 25L,
  cv_shift = 5L,
  intensity_scale = c(0.9, 1.1),
  noise_sd = 0,
  readout = c("gap", "attention", "presence"),
  presence_head = NULL,
  n_queries = 4L,
  seed = 1L,
  device = if (torch::cuda_is_available()) "cuda" else "cpu"
)
```

## Arguments

- Z_list:

  List of preprocessed, padded intensity matrices.

- y:

  Binary labels (integer or factor, two unique values).

- encoder:

  A pre-trained `dms_encoder`.

- k:

  Number of CV folds (default 5).

- epochs:

  Number of training epochs per fold.

- batch_size:

  Mini-batch size.

- lr:

  Adam learning rate.

- weight_decay:

  Adam L2 weight decay.

- dropout_p:

  Dropout probability before the linear classifier.

- freeze_backbone:

  If TRUE, encoder weights are not updated during fine-tuning.

- unfreeze_last_block:

  If TRUE, the last encoder residual block (`layer3`) is set trainable
  even when `freeze_backbone = TRUE`.

- augment:

  If TRUE, apply training-time augmentation
  ([`augment_batch()`](https://mljaniczek.github.io/classydms/reference/augment_batch.md):
  random RT/CV shifts, intensity scaling, noise).

- rt_shift, cv_shift:

  Maximum augmentation shift magnitudes (pixels) along RT and CV.

- intensity_scale:

  Length-2 numeric vector giving the bounds of random intensity
  rescaling.

- noise_sd:

  SD of Gaussian noise added to augmented images (0 disables).

- readout, presence_head, n_queries:

  Readout replacing GAP; see
  [`build_pretrained_classifier()`](https://mljaniczek.github.io/classydms/reference/build_pretrained_classifier.md).

- seed:

  RNG seed for fold assignment and training. Folds are drawn with this
  seed after the tensor dataset is built, so they do not depend on how
  much of R's RNG stream torch consumed on its lazy first-use
  initialisation. Fold `j` trains with seed `seed + 99 + j` (which is
  `100 + j` at the default `seed = 1`).

- device:

  Torch device for fine-tuning, e.g. `"cpu"`, `"cuda"`, or `"mps"`.
  Defaults to CUDA when available, else CPU (the previous behaviour).

## Value

List with `probs` (out-of-fold predictions) and `auc`.
