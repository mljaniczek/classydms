# Load a pre-trained encoder by name

Downloads the requested encoder file from GitHub Releases (cached
locally) and returns the loaded torch module. Pass
`with_autoencoder = TRUE` to also fetch the full autoencoder (encoder +
decoder), which is needed for reconstruction QC but not for
classification or saliency mapping.

## Usage

``` r
load_pretrained_encoder(
  name,
  repo = "mljaniczek/classydms",
  tag = "pretrained-v0",
  with_autoencoder = FALSE,
  force_download = FALSE
)
```

## Arguments

- name:

  Base name of the encoder file (without `.pt` extension).

- repo:

  GitHub repo.

- tag:

  Release tag.

- with_autoencoder:

  If TRUE, also fetch the autoencoder file and return both as a list.

- force_download:

  Force re-download even if cached locally.

## Value

Either a torch encoder module (if `with_autoencoder = FALSE`) or a list
with `encoder`, `autoencoder`, and `manifest` elements.

## Examples

``` r
if (FALSE) { # \dontrun{
enc <- load_pretrained_encoder("flu_v1")
bundle <- load_pretrained_encoder("flu_v1", with_autoencoder = TRUE)
bundle$encoder; bundle$autoencoder; bundle$manifest
} # }
```
