# Compound-anchored attention pooling

One learned query per catalog compound. Attention logits are
`q_k . key_t / sqrt(d) + log_prior_k(t)`, where the log-prior is a
Gaussian centred on compound k's catalog location with width taken from
[`compound_head_geometry()`](https://mljaniczek.github.io/classydms/reference/compound_head_geometry.md)
(floored at one feature-map cell, and scaled by a single learnable
factor). The prior tells each query WHERE to look; the content term
picks the peak within that window. Values are the raw encoder features,
so each compound's output token lives in the encoder's channel space.

## Usage

``` r
compound_attention_pool(in_channels, geometry, d_attn = 32L)
```

## Arguments

- in_channels:

  Encoder output channels.

- geometry:

  Output of
  [`compound_head_geometry()`](https://mljaniczek.github.io/classydms/reference/compound_head_geometry.md).

- d_attn:

  Query/key dimension.

## Value

Module whose forward returns
`list(tokens = (B, K, C), attn = (B, K, h*w))`.
