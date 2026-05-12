# Attach a classifier head to a pre-trained encoder

Deep-copies the encoder so each fold gets its own independent weights
(prevents cross-fold leakage when `unfreeze_last_block = TRUE`), wraps
it with adaptive average pooling and a dropout + linear classifier head.

## Usage

``` r
build_pretrained_classifier(
  encoder,
  num_classes = 2L,
  dropout_p = 0.3,
  freeze_backbone = TRUE,
  unfreeze_last_block = TRUE
)
```

## Arguments

- encoder:

  A pre-trained `dms_encoder`.

- num_classes:

  Number of output classes (default 2).

- dropout_p:

  Dropout probability before the linear layer.

- freeze_backbone:

  If TRUE, encoder weights are not updated.

- unfreeze_last_block:

  If TRUE, the encoder's last residual block (layer3) is set trainable
  even if `freeze_backbone = TRUE`.

## Value

An instantiated torch module.
