# Pre-training in depth

This article walks through the synthetic pre-training stage in detail:
estimating peak parameters from real data, generating synthetic samples
that match those distributions, training a denoising autoencoder via
self-supervised learning, and reading the resulting diagnostics. By the
end you should be able to choose `pretrain_*` arguments deliberately
rather than by trial and error.

This article assumes you have already followed
[`vignette("getting-started", package = "classydms")`](https://mljaniczek.github.io/classydms/articles/getting-started.md)
to load and preprocess your data into a list of padded `Z` matrices.

## Why pre-train?

Clinical GC-DMS cohorts typically contain 100–200 subjects. Training a
deep convolutional encoder on so few labeled examples leads to severe
overfitting: the network memorizes its training set rather than learning
generalizable features. Transfer learning sidesteps this by first
training the encoder on a large auxiliary dataset (where labels are not
needed), so that the resulting weights are a good starting point for the
actual classification task.

`classydms` does this with **synthetic data matched to your cohort**.
The synthetic generator builds 2-D Gaussian “peaks” with locations,
widths, and intensities sampled from distributions estimated empirically
from your real spectra. Pretraining the encoder to reconstruct or
denoise these synthetic images forces it to learn the morphology of
GC-DMS data — peaks against background, sparse signal in a large image,
characteristic spatial scales — without ever seeing a class label.

There are two pretraining entry points:

| Function | Memory | Total samples | When to use |
|----|----|----|----|
| [`pretrain_autoencoder()`](https://mljaniczek.github.io/classydms/reference/pretrain_autoencoder.md) (alias `pretrain_denoising`) | High — entire synthetic dataset in RAM | up to ~10,000 | Quick experiments, small synthetic datasets |
| [`pretrain_denoising_online()`](https://mljaniczek.github.io/classydms/reference/pretrain_denoising_online.md) | Constant — one batch at a time | hundreds of thousands to millions | Production runs, scaling up |

For real classification work the online variant is strongly recommended.
The rest of this article focuses on it.

## Step 1 — Estimate peak parameter distributions

``` r

library(classydms)
# Z_raw and Z_pretrim should be lists of trimmed sample matrices,
# parallel to your padded Z_padded list (see getting-started vignette)
peak_params <- estimate_peak_params(
  Z_raw,
  Z_pretrim_list = Z_pretrim,
  top_k = 150,
  frac_height = 0.5,    # FWHM-style width measurement
  verbose = TRUE
)
```

**What this does:** for each real sample, detect up to `top_k` local
maxima above the dust threshold, measure each peak’s width along RT and
CV at `frac_height` of peak height, and convert FWHM to Gaussian sigma.
Pool across all samples (cases and controls together — no labels used),
then fit log-normal distributions to the sigma and intensity vectors.
Also estimate background noise statistics from non-peak regions.

**Key fields in the return value:**

- `peak_params$n_peaks$mean` — average peaks detected per sample
  (synthetic images will draw their peak count near this value)
- `peak_params$sigma_rt$meanlog`, `sigma_rt$sdlog` — log-normal
  parameters of peak width along RT
- `peak_params$sigma_cv$meanlog`, `sigma_cv$sdlog` — same along CV
- `peak_params$intensity$meanlog`, `intensity$sdlog` — log-normal
  parameters of peak amplitude
- `peak_params$noise$mean`, `noise$sd` — additive background noise
  statistics
- `peak_params$rt_loc_raw`, `cv_loc_raw`, `sigma_rt_raw`,
  `sigma_cv_raw`, `intensity_raw` — raw observation vectors for
  diagnostic plotting

**Diagnostic check:** the verbose output should look something like

    sigma_rt range: 0.425 - 3.82 (median 1.91)
    sigma_cv range: 0.849 - 2.12 (median 1.27)

If the medians are far from one or two pixels, your peaks are either
much wider or much narrower than typical GC-DMS data and you may need to
revisit dust thresholding or trimming parameters.

## Step 2 — Visualize synthetic vs real samples (sanity check)

Before launching a long pretraining run, generate a few synthetic
samples and compare them to real ones by eye. This is the single most
important sanity check.

``` r

set.seed(7)
synth_examples <- lapply(1:3, function(i) {
  s <- generate_one_synthetic(peak_params,
                               H = pad_dims$H, W = pad_dims$W,
                               size_jitter = 0.6)
  normalize_sample(s$noisy)
})

# Side-by-side: top row = 3 real, bottom row = 3 synthetic
op <- par(mfrow = c(2, 3), mar = c(1, 1, 2, 1))
for (i in 1:3) {
  image(t(Z_padded[[i]])[, nrow(Z_padded[[i]]):1],
        col = hcl.colors(100, "YlOrRd", rev = TRUE),
        axes = FALSE, main = paste("Real #", i))
}
for (i in 1:3) {
  image(t(synth_examples[[i]])[, nrow(synth_examples[[i]]):1],
        col = hcl.colors(100, "YlOrRd", rev = TRUE),
        axes = FALSE, main = paste("Synth #", i))
}
par(op)
```

You want the two rows to look qualitatively similar in: number of peaks,
range of peak sizes, spatial distribution across RT and CV, and noise
level. If synthetic samples have far fewer peaks, peaks that are all the
same size, or peaks bunched in one region, revisit `peak_params` —
usually with `top_k`, `frac_height`, or by inspecting
`peak_params$rt_loc_raw` distributions.

## Step 3 — Run the online pretraining

``` r

pretrain_result <- pretrain_denoising_online(
  peak_params,
  H = pad_dims$H, W = pad_dims$W,
  steps_per_epoch  = 1000L,    # batches per epoch
  epochs           = 50L,
  batch_size       = 32L,
  lr               = 1e-3,
  weight_decay     = 1e-4,
  add_noise        = TRUE,     # denoising AE; FALSE = reconstruction AE
  size_jitter      = 0.6,
  grad_clip        = 1.0,
  norm_clamp       = 10.0,
  val_n            = 200L,     # fixed validation set size
  checkpoint_every = 10L,
  save_path        = "encoder.pt",
  seed             = 42L
)
```

### Choosing `steps_per_epoch` and `epochs`

Total samples seen by the network =
`steps_per_epoch * epochs * batch_size`. Targeting around **500,000 to
1,000,000 total samples** is a reasonable range for a 1,500 × 100 pixel
input space at default architecture. Anything below ~100,000 may not
fully exercise the encoder; above ~2,000,000 usually shows diminishing
returns.

Equivalent samples ↔︎ time trade-off:

| Goal             | steps_per_epoch | epochs | total |
|------------------|-----------------|--------|-------|
| Quick smoke test | 50              | 5      | 8 k   |
| Reasonable run   | 500             | 30     | 480 k |
| Production run   | 1000            | 50     | 1.6 M |
| Aggressive       | 2000            | 50     | 3.2 M |

### Outputs produced by the run

The function writes four families of files to disk when `save_path` is
set:

1.  `encoder.pt` — final encoder weights
2.  `encoder_autoencoder.pt` — full autoencoder (encoder + decoder),
    needed for reconstruction QC
3.  `encoder_manifest.Rdata` — bundles `peak_params`, all
    hyperparameters, loss history, validation loss history, batch-loss
    statistics, total samples seen, and timestamp. **Saved every
    epoch**, so a crash leaves usable artifacts.
4.  `encoder_epoch10.pt`, `encoder_epoch20.pt`, … — intermediate
    checkpoints (deletable after the run succeeds)

### Reading the per-epoch log

    Epoch 5/50 | train MSE: 0.0041 | val MSE: 0.0048
                | median: 0.0027 | p95: 0.008 | max: 0.078
                | Time: 29.0 min | ETA: 1305 min

- **train MSE** — average loss across all `steps_per_epoch` batches in
  this epoch. Each batch is fresh synthetic data, so this varies with
  the random sample composition.
- **val MSE** — loss on a fixed `val_n` synthetic samples held constant
  across epochs. This is the metric that is directly comparable across
  epochs and is what you should plot to assess convergence.
- **median / p95 / max** — distribution of *per-batch* losses within
  this epoch. Healthy: max within 5–10× median. Unhealthy: max ≫ median,
  indicating outlier batches the network is struggling on (re-check
  `peak_params` and `norm_clamp`).
- **ETA** — projected remaining time based on average per-epoch
  wall-clock so far.

### Reading the loss curve

``` r

load("encoder_manifest.Rdata")  # makes 'training_manifest' available

train_loss <- training_manifest$loss_history
val_loss   <- training_manifest$val_loss_history
epochs_x   <- seq_along(train_loss)

ylim <- range(c(train_loss, val_loss[is.finite(val_loss)]), na.rm = TRUE)
plot(epochs_x, train_loss, type = "l", lwd = 2,
     xlab = "Epoch", ylab = "MSE Loss", ylim = ylim,
     main = "Pre-Training Loss")
if (any(is.finite(val_loss))) {
  lines(epochs_x, val_loss, lwd = 2, col = "red")
  legend("topright",
         legend = c("Train (random batches)", "Val (fixed 200 samples)"),
         col = c("black", "red"), lwd = 2, bty = "n")
}
```

Three things to look for:

1.  **Both curves trending down together.** Healthy.
2.  **Train decreasing but val flat or increasing.** Overfitting — too
    many epochs for the current synthetic data quality. Stop earlier or
    scale up `steps_per_epoch`.
3.  **Both curves flat from the start.** Network can’t learn. Most
    common cause: pathological `peak_params` (e.g., almost zero peaks
    per sample, or extreme synthetic intensity range). Inspect Step 2
    sanity-check plots again.

### Why a fixed validation set matters

In online pretraining every training batch is freshly generated, so the
train loss is a moving target — comparing epoch 5’s train MSE to epoch
50’s train MSE is comparing losses on different synthetic samples. The
fixed validation set, generated once at the start of training with a
distinct seed and then held constant, gives a stable reference: a drop
in val MSE from epoch 5 to epoch 50 unambiguously means the encoder has
improved.

## Step 4 — Reload and use the encoder

The encoder is what transfers to downstream classification; the decoder
is discarded except for reconstruction QC.

``` r

encoder <- torch::torch_load("encoder.pt")
# Or, for reconstruction QC:
autoencoder <- torch::torch_load("encoder_autoencoder.pt")
# And the run metadata:
load("encoder_manifest.Rdata")
training_manifest$hyperparams      # how it was trained
training_manifest$peak_params      # the synthetic distribution
```

If you want to skip training entirely and use an encoder a collaborator
has shared via GitHub Releases:

``` r

encoder <- load_pretrained_encoder("flu_v1")
# Or with manifest + autoencoder for full reproducibility / QC:
bundle <- load_pretrained_encoder("flu_v1", with_autoencoder = TRUE)
```

## Step 5 — Reconstruction QC on real data

Feed real images through the autoencoder and inspect how well they’re
reconstructed. If reconstruction is good, the encoder has learned
features that transfer from synthetic to real. If reconstruction is
poor, the encoder has overfit to synthetic-specific quirks.

``` r

recon <- reconstruct_real_samples(autoencoder, Z_padded,
                                   sample_indices = 1:6)
op <- par(mfrow = c(2, 6), mar = c(1, 1, 2, 1))
for (r in recon) image(t(r$original)[, nrow(r$original):1],
                       col = hcl.colors(100, "YlOrRd", rev = TRUE),
                       axes = FALSE, main = paste0("Real #", r$index))
for (r in recon) image(t(r$reconstructed)[, nrow(r$reconstructed):1],
                       col = hcl.colors(100, "YlOrRd", rev = TRUE),
                       axes = FALSE,
                       main = sprintf("Recon (MSE=%.4f)", r$mse))
par(op)
```

Reconstructions should preserve the *positions* of major peaks even if
the *intensities* are softer (the bottleneck is lossy by design). Peaks
vanishing entirely is a warning sign. Reconstructions that look like a
smeared blob with no peak structure is a strong warning sign.

## Step 6 — Continue to downstream classification

Once you have a usable encoder, proceed to
[`vignette("interpretability", package = "classydms")`](https://mljaniczek.github.io/classydms/articles/interpretability.md)
for the classification + saliency mapping pipeline, or call
[`cv_encoder_elastic_net()`](https://mljaniczek.github.io/classydms/reference/cv_encoder_elastic_net.md)
directly for a quick AUC estimate.

## Resuming after a crash

Long training runs sometimes die mid-way (memory pressure, accidental
session restart, laptop sleep). The function writes the encoder,
autoencoder, and manifest after **every** epoch, so even an unexpected
interruption leaves usable artifacts on disk. To pick up where you left
off, pass the encoder path to `resume_from`:

``` r

pretrain_result <- pretrain_denoising_online(
  peak_params,
  H = pad_dims$H, W = pad_dims$W,
  steps_per_epoch = 1000L,
  epochs           = 50L,        # same as the original run
  batch_size       = 32L,
  save_path        = "encoder.pt",
  resume_from      = "encoder.pt",   # NEW: pick up from disk
  seed             = 42L
)
```

The function looks for `encoder_autoencoder.pt` and
`encoder_manifest.Rdata` next to the encoder file, loads the model
state, reads `last_epoch_completed` from the manifest, and continues
from the next epoch with the prior loss history preserved. If the
manifest says the run is already complete, the function returns
immediately without retraining.

**Caveat:** Adam optimizer moment estimates are not restored from disk,
so resumed training effectively re-warms Adam from scratch for a few
steps. For long runs this is invisible noise. The other hyperparameters
(`H`, `W`, `epochs`, `steps_per_epoch`, etc.) must match the original
run — change them and the resume becomes nonsense.

## Common pitfalls

- **Loss spikes to millions early in training.** The `norm_clamp` and
  `grad_clip` safeguards should prevent this. If you see it anyway, your
  `peak_params` likely produce pathologically extreme synthetic images.
  Inspect `peak_params$intensity$meanlog` and `peak_params$noise$sd` for
  sane values.
- **All encoder channels become identical at evaluation.** Mode collapse
  — the network gave up on most of its capacity. Often fixed by larger
  `steps_per_epoch`, more `size_jitter` to make synthetic peaks more
  diverse, or by adding noise (`add_noise = TRUE`).
- **`val MSE` reads `NA`.** The fixed validation set wasn’t built (this
  happens when `val_n = 0` is passed). Set `val_n > 0`.
