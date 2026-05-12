# List available pre-trained encoder files

Queries the GitHub Releases of `mljaniczek/classydms` and returns the
list of `.pt` and `.Rdata` files attached as release assets.

## Usage

``` r
list_pretrained_encoders(repo = "mljaniczek/classydms", tag = "pretrained-v0")
```

## Arguments

- repo:

  GitHub repo in `"owner/name"` form (default "mljaniczek/classydms").

- tag:

  Release tag to query (default "pretrained-v0"); use NULL to list all
  releases.

## Value

Character vector of available file names.
