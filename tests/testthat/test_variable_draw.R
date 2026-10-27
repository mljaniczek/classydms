# Variable-layer (and anchor) draw modes in generate_noisy_pretrain_sample().

# Reference captured on main (6d742aa) before the draw modes were added,
# then regenerated for the exclude_anchors fix: draw 2 drops anchor 2
# (prevalence 0.9, inside the variable range), which used to be re-drawn
# as a variable compound. Draws 1, 3 and 4 are identical to the 6d742aa
# capture. Regenerate only if uniform-mode output is meant to change.
uniform_reference_draws <- function() {
  cat_ <- make_toy_catalog()
  set.seed(123)
  lapply(1:4, function(i) {
    r <- generate_noisy_pretrain_sample(
      cat_, H = 64, W = 16, anchor_ids = c(1L, 2L),
      variable_config = list(n_per_sample = c(2, 5),
                             prevalence_range = c(0.05, 0.9)),
      contamination_config = list(n_per_sample = 3), anchor_dropout = 0.2)
    r[c("clean", "noisy", "anchor_ids", "variable_ids", "rendered_ids",
        "contamination_positions")]
  })
}

test_that("uniform draw (default and explicit) reproduces the pre-change output", {
  ref <- readRDS(test_path("fixtures", "uniform_reference.rds"))
  expect_identical(uniform_reference_draws(), ref)

  cat_ <- make_toy_catalog()
  run <- function(vc) {
    set.seed(7)
    generate_noisy_pretrain_sample(cat_, H = 64, W = 16, anchor_ids = 1L,
                                   variable_config = vc)
  }
  base <- list(n_per_sample = 3, prevalence_range = c(0, 1))
  expect_identical(run(base), run(c(base, draw = "uniform",
                                    prevalence_scale = 5)))
})

# Count how often each compound lands in the rendered variable layer.
variable_rates <- function(cat_, vc, n, anchor_ids = c(1L, 2L)) {
  ids <- cat_$compounds$compound_id
  hits <- matrix(FALSE, n, length(ids))
  labels_match <- TRUE
  for (k in seq_len(n)) {
    r <- generate_noisy_pretrain_sample(cat_, H = 64, W = 16,
                                        anchor_ids = anchor_ids,
                                        variable_config = vc,
                                        add_noise = FALSE)
    labels_match <- labels_match &&
      setequal(r$rendered_ids$variable, r$variable_ids)
    hits[k, ] <- ids %in% r$rendered_ids$variable
  }
  # Presence-head labels (rendered_ids) are exactly the drawn compounds.
  expect_true(labels_match)
  colnames(hits) <- ids
  hits
}

test_that("prevalence draw fires each compound at min(1, prevalence * scale)", {
  cat_ <- make_toy_catalog()
  prev <- cat_$compounds$prevalence
  anchors <- c(1L, 2L)
  n <- 1500
  for (scale in c(1, 1.5)) {
    set.seed(11)
    vc <- list(draw = "prevalence", prevalence_scale = scale,
               prevalence_range = c(0, 1))
    hits <- variable_rates(cat_, vc, n, anchor_ids = anchors)
    p <- pmin(1, prev * scale)
    eligible <- !(cat_$compounds$compound_id %in% anchors)
    rate <- colMeans(hits)
    # Anchors are excluded from the pool (exclude_anchors defaults TRUE).
    expect_true(all(rate[!eligible] == 0))
    # 4.5 binomial SDs: a false failure is very unlikely.
    tol <- 4.5 * sqrt(p * (1 - p) / n) + 1e-9
    expect_true(all(abs(rate[eligible] - p[eligible]) <= tol[eligible]),
                info = paste("scale", scale))
    pe <- p[eligible]
    expect_lte(abs(mean(rowSums(hits)) - sum(pe)),
               4.5 * sqrt(sum(pe * (1 - pe)) / n))
  }
})

test_that("prevalence draw honours prevalence_range, exclude_anchors and the cap", {
  cat_ <- make_toy_catalog()
  set.seed(3)
  hits <- variable_rates(cat_, list(draw = "prevalence",
                                    prevalence_range = c(0, 1),
                                    exclude_anchors = FALSE), 300)
  # Anchors 1, 2 (prevalence 0.95, 0.9) can now be drawn.
  expect_gt(mean(hits[, "1"]), 0.8)

  hits <- variable_rates(cat_, list(draw = "prevalence",
                                    prevalence_range = c(0.2, 0.7)), 300)
  out <- cat_$compounds$prevalence < 0.2 | cat_$compounds$prevalence > 0.7
  expect_true(all(colMeans(hits)[out] == 0))

  hits <- variable_rates(cat_, list(draw = "prevalence", prevalence_scale = 10,
                                    prevalence_range = c(0, 1),
                                    n_per_sample = 2), 50)
  expect_true(all(rowSums(hits) == 2))
})

test_that("prevalence draw keeps intensity_scale per layer", {
  cat_ <- make_toy_catalog()
  set.seed(5)
  r1 <- generate_noisy_pretrain_sample(
    cat_, H = 64, W = 16, anchor_ids = integer(0), add_noise = FALSE,
    location_jitter_rt = 0, location_jitter_cv = 0, size_jitter = 0,
    variable_config = list(draw = "prevalence", prevalence_scale = 100,
                           prevalence_range = c(0, 1),
                           intensity_scale = c(1, 1)))
  set.seed(5)
  r2 <- generate_noisy_pretrain_sample(
    cat_, H = 64, W = 16, anchor_ids = integer(0), add_noise = FALSE,
    location_jitter_rt = 0, location_jitter_cv = 0, size_jitter = 0,
    variable_config = list(draw = "prevalence", prevalence_scale = 100,
                           prevalence_range = c(0, 1),
                           intensity_scale = c(0.25, 0.25)))
  expect_length(r1$rendered_ids$variable, nrow(cat_$compounds))
  expect_equal(r2$layers$variable, 0.25 * r1$layers$variable)
})

test_that("anchor_draw = 'prevalence' fires anchors at prevalence * (1 - dropout)", {
  cat_ <- make_toy_catalog()
  anchors <- c(3L, 8L)                       # prevalence 0.6, 0.8
  p <- cat_$compounds$prevalence[anchors]
  n <- 1500
  for (drop in c(0, 0.25)) {
    set.seed(9)
    hits <- t(vapply(seq_len(n), function(k) {
      r <- generate_noisy_pretrain_sample(cat_, H = 64, W = 16,
                                          anchor_ids = anchors,
                                          anchor_draw = "prevalence",
                                          anchor_dropout = drop,
                                          add_noise = FALSE)
      stopifnot(setequal(r$rendered_ids$anchors, r$anchor_ids))
      anchors %in% r$rendered_ids$anchors
    }, logical(2)))
    q <- p * (1 - drop)
    expect_true(all(abs(colMeans(hits) - q) <= 4.5 * sqrt(q * (1 - q) / n)))
  }
})

# Configs the exclude_anchors fix must not change: no thinning, or
# thinning with a variable range below every anchor's prevalence.
# Reference captured on main (df9a782) before the fix.
unchanged_case_draws <- function() {
  cat_ <- make_toy_catalog()
  keep <- c("clean", "noisy", "anchor_ids", "variable_ids", "rendered_ids")
  one <- function(seed, ...) {
    set.seed(seed)
    lapply(1:3, function(i)
      generate_noisy_pretrain_sample(cat_, H = 64, W = 16,
                                     anchor_ids = c(1L, 2L, 8L), ...)[keep])
  }
  list(
    none_uniform = one(1, variable_config = list(n_per_sample = c(2, 5),
                                                 prevalence_range = c(0.05, 1))),
    none_prev    = one(2, variable_config = list(draw = "prevalence",
                                                 prevalence_range = c(0.05, 1))),
    below_drop   = one(3, anchor_dropout = 0.4,
                       variable_config = list(n_per_sample = c(2, 5),
                                              prevalence_range = c(0.05, 0.6))),
    below_prev   = one(4, anchor_draw = "prevalence",
                       variable_config = list(draw = "prevalence",
                                              prevalence_range = c(0.05, 0.6))))
}

test_that("exclude_anchors fix leaves unthinned / out-of-range configs unchanged", {
  ref <- readRDS(test_path("fixtures", "unchanged_anchor_exclusion.rds"))
  expect_identical(unchanged_case_draws(), ref)
})

# Run n draws with anchors whose prevalence sits inside the variable
# range; return per-anchor firing indicators across both layers.
anchor_leak_rates <- function(cat_, anchors, n, vc, ...) {
  leaked <- 0L
  hits <- t(vapply(seq_len(n), function(k) {
    r <- generate_noisy_pretrain_sample(cat_, H = 64, W = 16,
                                        anchor_ids = anchors,
                                        variable_config = vc,
                                        add_noise = FALSE, ...)
    leaked <<- leaked + sum(r$variable_ids %in% anchors) +
      sum(r$rendered_ids$variable %in% anchors)
    stopifnot(setequal(r$rendered_ids$anchors, r$anchor_ids))
    anchors %in% unlist(r$rendered_ids)
  }, logical(length(anchors))))
  expect_equal(leaked, 0L)
  colMeans(hits)
}

test_that("anchors thinned by anchor_draw = 'prevalence' are not re-drawn", {
  cat_ <- make_toy_catalog()
  anchors <- c(1L, 2L, 8L)                   # prevalence 0.95, 0.9, 0.8
  p <- cat_$compounds$prevalence[anchors]
  n <- 2000
  set.seed(21)
  rate <- anchor_leak_rates(cat_, anchors, n, anchor_draw = "prevalence",
                            vc = list(draw = "prevalence",
                                      prevalence_range = c(0.05, 1),
                                      intensity_scale = c(1, 1)))
  # Before the fix these fired at 1 - (1 - p)^2 (0.9975, 0.99, 0.96).
  expect_true(all(abs(rate - p) <= 4.5 * sqrt(p * (1 - p) / n)))
})

test_that("anchors silenced by anchor_dropout are not re-drawn (uniform)", {
  cat_ <- make_toy_catalog()
  anchors <- c(1L, 2L, 8L)
  n <- 2000
  q <- rep(1 - 0.3, length(anchors))
  set.seed(22)
  rate <- anchor_leak_rates(cat_, anchors, n, anchor_dropout = 0.3,
                            vc = list(n_per_sample = c(3, 6),
                                      prevalence_range = c(0.05, 1)))
  # Before the fix these fired at ~0.855.
  expect_true(all(abs(rate - q) <= 4.5 * sqrt(q * (1 - q) / n)))
})

test_that("bad draw settings fail clearly", {
  cat_ <- make_toy_catalog()
  gen <- function(vc, ...) generate_noisy_pretrain_sample(
    cat_, H = 64, W = 16, anchor_ids = 1L, variable_config = vc, ...)
  expect_error(gen(list(draw = "weighted")), "must be \"uniform\" or \"prevalence\"")
  expect_error(gen(list(draw = c("uniform", "prevalence"))), "must be")
  expect_error(gen(list(draw = "prevalence", prevalence_scale = -1)),
               "prevalence_scale must be")
  expect_error(gen(list(draw = "prevalence", prevalence_scale = NA_real_)),
               "prevalence_scale must be")
  expect_error(gen(NULL, anchor_draw = "sometimes"))

  no_prev <- cat_
  no_prev$compounds$prevalence <- NULL
  expect_error(generate_noisy_pretrain_sample(
    no_prev, H = 64, W = 16, anchor_ids = 1L,
    variable_config = list(draw = "prevalence")), "needs a `prevalence`")
  expect_error(generate_noisy_pretrain_sample(
    no_prev, H = 64, W = 16, anchor_ids = 1L, anchor_draw = "prevalence"),
    "needs a `prevalence`")
})

test_that("pretrain_denoising_online validates the draw mode at entry", {
  cat_ <- make_toy_catalog()
  expect_error(pretrain_denoising_online(
    catalog = cat_, H = 64, W = 16, anchor_ids = 1L,
    variable_config = list(draw = "weighted"), verbose = FALSE),
    "must be \"uniform\" or \"prevalence\"")
  expect_error(pretrain_denoising_online(
    catalog = cat_, H = 64, W = 16, anchor_ids = 1L,
    variable_config = list(draw = "prevalence", prevalence_scale = -2),
    verbose = FALSE), "prevalence_scale must be")
})

test_that("prevalence draw is valid and independent across parallel workers", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  cat_ <- make_toy_catalog()
  vc <- list(draw = "prevalence", prevalence_range = c(0, 1))
  old_plan <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old_plan), add = TRUE)
  # Workers are fresh R sessions that load the *installed* classydms
  # (under devtools::test() that may be an older build).
  worker_args <- future::value(future::future(tryCatch(
    names(formals(classydms::generate_noisy_pretrain_sample)),
    error = function(e) character(0))))
  skip_if_not("anchor_draw" %in% worker_args,
              "installed classydms predates the draw modes")
  run <- function() future.apply::future_lapply(seq_len(40), function(i)
    classydms::generate_noisy_pretrain_sample(cat_, H = 64, W = 16,
                                              anchor_ids = c(1L, 2L),
                                              variable_config = vc),
    future.seed = 2024L)
  a <- run(); b <- run()
  expect_identical(a, b)                     # reproducible under a seed
  ids <- vapply(a, function(r) paste(sort(r$rendered_ids$variable),
                                     collapse = ","), "")
  expect_gt(length(unique(ids)), 10)         # streams differ per task
  expect_true(all(vapply(a, function(r)
    all(is.finite(r$noisy)) && !any(r$rendered_ids$variable %in% 1:2), NA)))
  p <- pmin(1, cat_$compounds$prevalence[-(1:2)])
  expect_equal(mean(lengths(lapply(a, function(r) r$rendered_ids$variable))),
               sum(p), tolerance = 0.25)
})

test_that("pretrain_denoising_online records the draw modes in the manifest", {
  skip_on_cran()
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed())
  cat_ <- make_toy_catalog()
  dir <- tempfile("pretrain_draw_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  sp <- file.path(dir, "enc.pt")
  suppressMessages(pretrain_denoising_online(
    catalog = cat_, H = 64, W = 32, anchor_ids = c(1L, 2L),
    variable_config = list(draw = "prevalence", prevalence_scale = 0.5,
                           prevalence_range = c(0, 1)),
    anchor_draw = "prevalence", presence_head = TRUE,
    steps_per_epoch = 1L, epochs = 1L, batch_size = 2L, val_n = 2L,
    save_path = sp, device = "cpu", verbose = FALSE))
  e <- new.env()
  load(sub("\\.pt$", "_manifest.Rdata", sp), envir = e)
  hp <- e$training_manifest$hyperparams
  expect_identical(hp$anchor_draw, "prevalence")
  expect_identical(hp$variable_draw, "prevalence")
  expect_identical(hp$variable_prevalence_scale, 0.5)
})
