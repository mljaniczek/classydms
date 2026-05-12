# Getting started with classydms

This vignette walks through the full `classydms` pipeline at a high
level: load raw GC-DMS data, preprocess, generate synthetic data,
pretrain the encoder, classify with elastic net, and produce a spatial
saliency map for biomarker localization.

``` r

library(classydms)
library(purrr)   # for map(), used in some downstream steps
```

## 1. Load raw data from a directory

The expected file layout is one file per subject, with disease and
control samples in separate folders. Each file should be tab-delimited
with the GC-DMS layout described in
[`?read_dms_file`](https://mljaniczek.github.io/classydms/reference/read_dms_file.md)
(line 1 = “Vc”, line 2 = tab + CV values, line 3 = header, lines 4+ =
time + intensities).

``` r

# One-shot loader: reads all files matching the POS-channel pattern,
# skips broken ones, names each sample after its file basename (so
# results stay traceable to subjects), and returns a (samples, y) pair.
disease_dir = here::here("sampledata/influenza")
control_dir = here::here("sampledata/controls")

loaded <- load_dms_directory(
  disease_dir = disease_dir,
  control_dir = control_dir,
  pattern     = "POS.*\\.(xls|txt|tsv|csv)$"  # POS channel only — see note
)

samples <- loaded$samples   # named list, each element has $path, $time, $cv, $Z
y       <- loaded$y         # named integer vector: 2 = disease, 1 = control

# Sample names map to file basenames without extensions, so you can:
samples[["241115-S1994-S1-DMS1_POS"]]    # index by subject
names(samples)                            # see the full list of subjects

# Sanity check
stopifnot(is.list(samples), length(samples) > 0,
          is.matrix(samples[[1]]$Z),
          all(nzchar(names(samples))))   # no empty-string names
```

**Channel choice**: classifiers in this package train on the **positive
ion channel only** (POS files). Negative-channel data is not used. The
default `pattern` enforces this — pass an explicit `pattern` argument
only if you need a different channel or filename convention.

If you have a non-standard file format, write your own loader that
produces the same list shape — each element must have `$path`, `$time`,
`$cv`, and `$Z`, and the outer list should have meaningful (non-empty)
names. As long as that contract is met, the rest of the pipeline works
unchanged.

## 2. Preprocess

``` r

# Savitzky-Golay smoothing + AsLS baseline correction + clip below zero
samples_proc <- lapply(samples, process_one_sample,
                        sg_p = 3, sg_n = 21,
                        asls_lambda = 1e6, asls_p = 0.01, asls_niter = 10)

# Dust threshold each sample to a fixed floor (e.g. global value, or
# per-sample via robust_eps). For a fixed floor:
dust_thr <- 0.005
samples_dust <- lapply(samples_proc, baseline_basement,
                        basement_thr = dust_thr)

# Cohort-wide trim bounds
b_list <- lapply(samples_dust,
                  function(s) trim_bounds_from_occupancy(s$Z))
R0 <- max(sapply(samples_dust, function(s) nrow(s$Z)))
C0 <- max(sapply(samples_dust, function(s) ncol(s$Z)))
global_b <- pooled_trim_bounds(b_list, R = R0, C = C0,
                                q_lo = 0.05, q_hi = 0.95)
samples_trimmed <- lapply(samples_dust, trim_sample, bounds = global_b)

# Per-sample log-quantile normalization + pad to common size
Z_norm   <- map(samples_trimmed, ~normalize_sample(.x$Z))
pad_dims <- compute_pad_targets(Z_norm)
Z_padded <- map(Z_norm, ~pad_to_target(.x, pad_dims$H, pad_dims$W))
```

## 3. Estimate peak parameters from real data (label-agnostic)

``` r

Z_raw     <- map(samples_trimmed, ~.x$Z)
Z_pretrim <- map(samples_trimmed,
                  ~if (!is.null(.x$Z_pretrim)) .x$Z_pretrim else .x$Z)

peak_params <- estimate_peak_params(Z_raw,
                                     Z_pretrim_list = Z_pretrim,
                                     top_k = 150,
                                     frac_height = 0.5)
```

## 4. Pretrain the encoder

For most users, the online variant is recommended because it removes the
RAM ceiling and lets the model see hundreds of thousands of synthetic
samples:

``` r

pretrain_result <- pretrain_denoising_online(
  peak_params, H = pad_dims$H, W = pad_dims$W,
  steps_per_epoch = 10L, epochs = 10L, batch_size = 32L,
  save_path = "encoder.pt", seed = 42L)
```

This writes `encoder.pt`, `encoder_autoencoder.pt`, and
`encoder_manifest.Rdata` to the working directory, plus periodic
checkpoints. The manifest contains the full hyperparameter configuration
and loss curves for reproducibility.

Alternatively, skip pretraining entirely by downloading an encoder a
collaborator has already trained:

``` r

encoder <- load_pretrained_encoder("flu_v1")
```

## 5. Classify with elastic net on encoder features

``` r

elnet <- cv_encoder_elastic_net(pretrain_result$encoder,
                                 Z_list = Z_padded, y = y,
                                 k = 5, alpha = 0.5)
elnet$auc
elnet$coefs
```

## 6. Saliency maps

``` r

disease_idx <- which(y == 2)
orig_dims <- lapply(samples_trimmed,
                     function(s) c(H = nrow(s$Z), W = ncol(s$Z)))

sal_disease <- aggregate_saliency_masked(pretrain_result$encoder,
                                          Z_padded, coefs = elnet$coefs,
                                          sample_indices = disease_idx,
                                          orig_dims = orig_dims)
```

The result is a 2-D matrix in padded image coordinates. Map back to
physical (RT in seconds, CV in volts) using `samples_trimmed[[1]]$time`
and `samples_trimmed[[1]]$cv`. Saliency hotspots correspond to (RT, CV)
locations whose encoder activation, weighted by elastic net
coefficients, most strongly drives the disease prediction. These
coordinates are the candidate biomarker compounds.

## Troubleshooting

### “attempt to use zero-length variable name”

This error means a list somewhere in your pipeline has an
**empty-string** name (`""`), and a downstream operation is trying to
use that name. Two common causes:

1.  The `samples` list was built via
    [`purrr::map()`](https://purrr.tidyverse.org/reference/map.html) +
    [`purrr::keep()`](https://purrr.tidyverse.org/reference/keep.html)
    (or [`c()`](https://rdrr.io/r/base/c.html) of filtered lists) and
    inherited empty-string names. **Fix**: replace any empty names with
    meaningful ones, or drop them entirely.
    [`load_dms_directory()`](https://mljaniczek.github.io/classydms/reference/load_dms_directory.md)
    names samples after their file basenames and explicitly verifies no
    name is empty.
2.  A sample’s `Z` matrix has empty-string column or row names. **Fix**:
    `dimnames(Z) <- NULL`. The
    [`read_dms_file()`](https://mljaniczek.github.io/classydms/reference/read_dms_file.md)
    function does this for you.

If you wrote a custom loader, give each sample a meaningful name
(e.g. file basename, subject ID) and make sure the `Z` matrices don’t
carry empty dimnames. Keeping non-empty names is important — you’ll want
them for tracing classification results and saliency maps back to
specific subjects.
