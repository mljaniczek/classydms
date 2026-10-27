# Attach a classifier head to a pre-trained encoder

Deep-copies the encoder so each fold gets its own independent weights
(prevents cross-fold leakage when `unfreeze_last_block = TRUE`), then
adds a readout that turns the `(B, C, h, w)` feature map into a vector,
followed by dropout + linear classifier.

## Usage

``` r
build_pretrained_classifier(
  encoder,
  num_classes = 2L,
  dropout_p = 0.3,
  freeze_backbone = TRUE,
  unfreeze_last_block = TRUE,
  readout = c("gap", "attention", "presence"),
  presence_head = NULL,
  n_queries = 4L
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

- readout:

  One of `"gap"` (default), `"attention"`, `"presence"`.

- presence_head:

  Required when `readout = "presence"`: the `presence_head` returned by
  [`pretrain_denoising_online()`](https://mljaniczek.github.io/classydms/reference/pretrain_denoising_online.md)
  (or loaded from its `*_presence_head.pt`).

- n_queries:

  Number of queries for `readout = "attention"`.

## Value

An instantiated torch module.

## Details

Readouts:

- `"gap"`: global average pooling (the original readout). Averages over
  all of RT x CV, so it discards where activation occurred.

- `"attention"`:
  [`attention_pool_readout()`](https://mljaniczek.github.io/classydms/reference/attention_pool_readout.md)
  — `n_queries` learned queries with 2D positional encoding, trained
  end-to-end. Works with any existing encoder.

- `"presence"`: a pretrained
  [`dms_presence_head()`](https://mljaniczek.github.io/classydms/reference/dms_presence_head.md)
  producing one logit per catalog compound; the classifier is then a
  linear model over compound presence evidence. The head is deep-copied
  and follows `freeze_backbone` (frozen by default).
