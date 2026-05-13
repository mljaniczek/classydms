# Compute, memory, and timing reference

This article gives expected wall-clock times, memory footprints, and
optimization knobs for
[`pretrain_denoising_online()`](https://mljaniczek.github.io/classydms/reference/pretrain_denoising_online.md)
across different hardware. Use it to plan training runs, choose
hyperparameters that fit your machine, and decide between CPU and GPU
backends.

The numbers below come from real runs on the influenza sample data with
padded dimensions `1472 x 96`. If your padded dimensions differ
substantially, scale linearly with `H * W`.

## The training budget

Total samples seen by the encoder =
`steps_per_epoch * epochs * batch_size`. For a useful encoder, plan for
somewhere between 500,000 and 2,000,000 total samples:

| Budget | Use case |
|----|----|
| 8,000 (50 steps × 5 epochs × 32) | Smoke test — verify pipeline works |
| 100,000 (500 × 10 × 32 ish) | Quick experiment to compare hyperparameters |
| 500,000 (1000 × 16 × 32 ish) | Reasonable production run |
| 1,000,000 (1000 × 32 × 32 ish) | Strong production run |
| 2,000,000 (1000 × 64 × 32 ish) | Aggressive — diminishing returns expected past here |

Going beyond 2M usually doesn’t help — the autoencoder has seen enough
variations to saturate. If it’s still improving past 1M, your
`peak_params` are probably under-fitted (too few peaks, too tight a
sigma distribution); fix that first rather than throwing more compute at
it.

## Backend matters a lot

Three backends are supported via the `device` argument:

- `"cpu"` — works everywhere; multi-threaded via libtorch’s internal
  pool
- `"mps"` — Metal Performance Shaders on Apple Silicon; **strongly
  recommended on M-series Macs**
- `"cuda"` — NVIDIA GPU; works if you have one and an appropriate
  libtorch installation

By default the function picks `"cuda"` if available, otherwise `"cpu"`.
On Apple Silicon you need to explicitly set `device = "mps"`. Verify MPS
is available with
[`torch::backends_mps_is_available()`](https://torch.mlverse.org/docs/reference/backends_mps_is_available.html).

The wall-clock difference between `"cpu"` and `"mps"` on the same
M-series chip is typically **10-15x in MPS’s favor**. Use MPS unless you
have a specific reason not to.

## Wall-clock estimates

All numbers assume padded `H = 1472, W = 96`, `batch_size = 32`, and
roughly default hyperparameters. Times are per-epoch and per
full-1M-sample run (~32 epochs × 1000 steps).

### Apple Silicon (MPS)

**Important caveat measured empirically**: on M-series GPUs, the
autoencoder forward/backward pass is so fast that the dominant
wall-clock cost is **R-side synthetic data generation**, not GPU
compute. The numbers below reflect this. If parallel data generation is
added in a future version (see `CLAUDE.md` “Potential optimizations”),
MPS times will drop by another 2-3x and the table below will need
updating.

| Chip | GPU cores | Per-batch time | 1 epoch (500 steps, batch 64) | 1M samples (32 epochs) |
|----|----|----|----|----|
| M1 Pro | 14-16 | ~6 sec | ~50 min | ~27 hours |
| M1 Max | 24-32 | ~5 sec | ~42 min | ~22 hours |
| M2 Pro | 16-19 | ~4.5 sec | ~38 min | ~20 hours |
| M2 Max | 30-38 | ~4 sec | ~33 min | ~18 hours |
| M3 Pro | 14-18 | ~4 sec | ~33 min | ~18 hours |
| M3 Max | 30-40 | ~3.7 sec | ~31 min | ~17 hours |
| M4 Pro | 16-20 | ~3.6 sec | ~30 min | ~16 hours |
| **M4 Max** | **32-40** | **~3.4 sec** | **~28 min** | **~15 hours** |

A 1M-sample run is an overnight job; a 2M-sample run is closer to 24-30
hours and benefits from being split across two nights via `resume_from`.
These numbers assume `batch_size = 64`. Larger batches help amortize the
per-batch R-side generation cost: doubling batch size from 32 to 64
roughly halves the per-sample time. `batch_size = 128` may give another
modest gain on M-Pro/Max chips; diminishing returns past that.

**Why the chip matters less than you’d expect**: the GPU portion of wall
time is 1-1.5 seconds per batch even on an M4 Max, while R-side data
generation is ~2 seconds. A faster GPU only speeds up the GPU half, so
going from M1 Pro to M4 Max saves ~1 second per batch, not ~10x.

### CPU only

For CPU runs, performance scales primarily with the number of
performance cores (E-cores contribute less per thread). Set
`num_threads` to roughly the number of performance cores on your
machine.

| Hardware                          | `num_threads` | 1 epoch | 1M samples |
|-----------------------------------|---------------|---------|------------|
| MacBook Air (M1, 4P+4E)           | 4             | ~75 min | ~40 hours  |
| MacBook Pro (M3 Pro, 6P+6E)       | 6             | ~45 min | ~24 hours  |
| MacBook Pro (M4 Pro, 10P+4E)      | 10            | ~30 min | ~16 hours  |
| MacBook Pro (M4 Max, 12P+4E)      | 12            | ~25 min | ~13 hours  |
| Linux workstation (8 cores, x86)  | 8             | ~50 min | ~27 hours  |
| Linux workstation (16 cores, x86) | 16            | ~30 min | ~16 hours  |

CPU runs are roughly **10-15x slower than MPS** on the same Apple
Silicon chip. If you have an M-series Mac, there’s almost no reason to
run on CPU.

### NVIDIA CUDA

| GPU          | VRAM     | 1 epoch  | 1M samples |
|--------------|----------|----------|------------|
| RTX 3060     | 12 GB    | ~3 min   | ~1.7 hours |
| RTX 3080     | 10-12 GB | ~2 min   | ~1.1 hours |
| RTX 4080     | 16 GB    | ~1.5 min | ~50 min    |
| RTX 4090     | 24 GB    | ~1 min   | ~35 min    |
| A100 (cloud) | 40-80 GB | ~45 sec  | ~25 min    |

CUDA is the fastest option overall. With plenty of VRAM, you can also
push `batch_size` to 128 or 256 for additional throughput.

## Memory footprint

The package is designed to run on machines with as little as 8 GB RAM,
but more is better for safety margin.

### Per-batch memory (clean + noisy tensors)

    mem_per_batch_MB = batch_size * 1 * H * W * 4 * 2 / 1e6

At default `batch_size = 32`, `H = 1472`, `W = 96`:

    32 * 1 * 1472 * 96 * 4 * 2 / 1e6 = 36 MB per batch

This is small. The bottleneck is not the batch tensors themselves but
the activations during the forward pass through the autoencoder.

### Activation memory during forward/backward

Roughly 5-10x the per-batch input size for a moderate-depth
encoder/decoder. So ~200-400 MB during the forward pass. With gradients
and Adam moments, total training memory is typically 1-2 GB peak.

### Validation set

    val_mem_MB = val_n * 1 * H * W * 4 * 2 / 1e6

At `val_n = 200` and default dims: ~225 MB. Held resident on the device
throughout training.

### Total expected memory during training

| Backend | Total peak memory                                    |
|---------|------------------------------------------------------|
| CPU     | 2-4 GB R session, plus libtorch internal allocations |
| MPS     | 2-4 GB unified memory (shared with system)           |
| CUDA    | 1-2 GB VRAM at batch_size 32, more at larger batches |

If you see memory pressure (yellow/red in Activity Monitor on macOS, or
`nvidia-smi` showing \>80% VRAM), drop `batch_size` or `val_n`.

## Optimization knobs in order of impact

1.  **Use MPS (or CUDA) instead of CPU** if available. By far the
    largest single speedup. 10-15x.

2.  **Increase `batch_size`** on a GPU backend if utilization is below
    80%. Bigger batches amortize kernel launch overhead. Halve
    `steps_per_epoch` to keep total samples constant. Diminishing
    returns past `batch_size = 128`.

3.  **Tune `num_threads`** on CPU runs. On Apple Silicon, set to the
    number of performance cores (typically 4-12 depending on chip). On
    x86 Linux, set to the number of physical cores (not threads —
    hyperthreading rarely helps for this workload).

4.  **Reduce `H * W`** if your data allows. Smaller padded dimensions
    linearly reduce per-batch compute. Achievable via more aggressive
    trimming or by accepting more padding by setting a smaller
    `compute_pad_targets(..., multiple = 16)`.

5.  **Drop `val_n` to 0**. Disables the fixed-validation-set eval each
    epoch. Saves a few seconds per epoch but loses the stable progress
    metric. Not recommended.

## Estimating wall-clock for your specific machine

Run a single-epoch smoke test before committing to a long run:

``` r

system.time({
  pretrain_result <- pretrain_denoising_online(
    peak_params, H = pad_dims$H, W = pad_dims$W,
    steps_per_epoch = 100L,
    epochs = 1L,
    batch_size = 32L,
    val_n = 50L,
    device = "mps",  # or "cpu" / "cuda"
    save_path = NULL,
    seed = 42L
  )
})
```

If this reports an elapsed time of T seconds for 100 steps, then for a
full 1M-sample run:

    estimated_hours = T * (1e6 / 32 / 100) / 3600
                    = T * 312.5 / 3600
                    = T * 0.087

A 100-step smoke test that takes 60 seconds implies a 1M run will take
~5.2 hours.

## What dominates wall time?

The relative contribution of each pipeline stage depends on the backend:

| Stage                                          | CPU    | MPS    | CUDA   |
|------------------------------------------------|--------|--------|--------|
| Synthetic data generation (R loop)             | 30-40% | 50-60% | 60-70% |
| Forward + backward pass                        | 50-60% | 30-40% | 25-35% |
| Optimizer step                                 | 5-10%  | 3-5%   | 2-3%   |
| Disk I/O (save_manifest, encoder, autoencoder) | 1-2%   | 1-2%   | 1-2%   |
| Validation forward pass                        | 1-2%   | 1-2%   | 1-2%   |

On CPU, the model forward/backward pass dominates; on GPU backends, the
**CPU-side synthetic data generation** becomes the bottleneck because
the model is fast. This is why parallel synthetic generation (see
“Potential optimizations” in `CLAUDE.md`) would help more on GPU than on
CPU.

## Activity Monitor / nvidia-smi checklist while training

### On macOS with MPS

Open Activity Monitor → Window → GPU History (`⌘+4`):

- `rsession-arm64` should appear consistently in the GPU users list at
  **40-80% GPU usage**.
- Memory tab should stay green throughout. Yellow or red = pressure,
  drop `batch_size`.
- Swap Used at the bottom of the Memory tab should stay near 0.
- CPU tab will show R using maybe 100-300% (i.e., 1-3 cores worth) for
  the synthetic data generation between batches. This is expected.

### On Linux with CUDA

Run `watch -n 1 nvidia-smi` in another terminal:

- GPU utilization should be **70-95%** during training.
- VRAM should be stable across epochs, not creeping up (that would
  indicate a memory leak).
- Power draw should be near the GPU’s rated TDP — if it’s at 50% TDP,
  your CPU or batch size is bottlenecking the GPU.

### Red flags

- **Memory pressure climbing across epochs**: indicates a leak. The
  [`gc()`](https://rdrr.io/r/base/gc.html) between epochs should prevent
  this, but if you see it, file an issue.
- **GPU utilization below 30%**: data generation is the bottleneck. Try
  increasing `batch_size` or accepting CPU-only execution.
- **Wall-clock per epoch increasing over time**: usually a memory issue
  causing swap or GC thrash. Reduce batch size or close other
  applications.

## Crash recovery

Long-running pretraining can fail mid-run for any number of reasons (OS
kill, MPS instability, RStudio session quirks, laptop sleep). The
function writes the encoder, autoencoder, and manifest **after every
epoch**, so the most recent completed epoch’s state is always on disk.

To resume from where it crashed:

``` r

pretrain_denoising_online(
  ...,                          # same hyperparameters as before
  save_path   = "encoder.pt",   # same path
  resume_from = "encoder.pt"    # NEW
)
```

The function reads the manifest, picks up at `last_epoch_completed + 1`,
and inherits prior loss history. If the manifest says the run is already
complete, it returns immediately without retraining.

Caveat: Adam optimizer moment estimates are NOT restored. The first few
resumed steps effectively re-warm Adam. For runs of 30+ epochs this is
invisible noise.

## When to choose what

A flowchart for picking hyperparameters and backend:

1.  **Apple Silicon Mac?** Use MPS. Set `device = "mps"`.
2.  **NVIDIA GPU with at least 8 GB VRAM?** Use CUDA. Set
    `device = "cuda"`. Push `batch_size` to 64-128.
3.  **Neither, but at least 8 CPU cores?** Use CPU. Set `num_threads` to
    the number of performance/physical cores.
4.  **A laptop with \< 8 GB RAM?** Drop `batch_size` to 16 or 8. Keep
    total samples modest (\< 500k).

For total sample budget:

1.  **First time running on a new dataset?** Start with 100k samples,
    sanity-check the val loss is decreasing.
2.  **Reasonable production run?** 500k-1M samples.
3.  **Best possible encoder?** 1M-2M samples.
4.  **Past 2M and still want more?** Improve `peak_params` quality
    first. Larger `top_k`, smaller `frac_height`, longer training-data
    tails — these typically help more than additional synthetic samples
    beyond 2M.
