# Compound-presence head

Compound-anchored attention pooling
([`compound_attention_pool()`](https://mljaniczek.github.io/classydms/reference/compound_attention_pool.md))
followed by a per-compound linear readout, producing one presence logit
per catalog compound. Trained as an auxiliary pretraining task in
[`pretrain_denoising_online()`](https://mljaniczek.github.io/classydms/reference/pretrain_denoising_online.md)
(`presence_head = TRUE`) against the generator's ground-truth list of
compounds it rendered. Downstream, the logits form an `(N, K)` feature
matrix whose columns ARE catalog compounds — see
[`extract_encoder_features()`](https://mljaniczek.github.io/classydms/reference/extract_encoder_features.md)
with `readout = "presence"`.

## Usage

``` r
dms_presence_head(in_channels, geometry, d_attn = 32L)
```

## Arguments

- in_channels:

  Encoder output channels.

- geometry:

  Output of
  [`compound_head_geometry()`](https://mljaniczek.github.io/classydms/reference/compound_head_geometry.md).

- d_attn:

  Query/key dimension.
