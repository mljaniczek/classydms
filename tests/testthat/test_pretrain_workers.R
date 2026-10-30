# Synthetic-data generation in pretrain_denoising_online(): serial and
# parallel paths (R/pretrain_workers.R).

# Batch digests captured on main at 9ab30cd, before generation moved to
# package-level worker functions (helper-pretrain-batches.R).
ref_batches <- function() readRDS(test_path("fixtures", "pretrain_batches_9ab30cd.rds"))

test_that("serial batches are unchanged from 9ab30cd", {
  skip_on_cran()
  expect_identical(batch_digests(capture_pretrain_batches(1L)),
                   ref_batches()$serial)
})

test_that("per-sample seeds match future_lapply(future.seed = TRUE)", {
  skip_if_not_installed("future.apply")
  old_plan <- future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
  for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
    RNGkind(kind)
    set.seed(3)
    ours <- .pretrain_rng_seeds(7L)
    ours_after <- .Random.seed
    set.seed(3)
    theirs <- future.apply::future_lapply(1:7, function(i) .Random.seed,
                                          future.seed = TRUE)
    expect_identical(ours, theirs)
    expect_identical(ours_after, .Random.seed)   # same main-stream advance
    expect_identical(RNGkind()[1], kind)
  }
})

test_that("parallel batches are reproducible, distinct and unchanged from 9ab30cd", {
  skip_if_workers_differ()
  old_plan <- future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)
  a <- capture_pretrain_batches(2L)
  expect_identical(capture_pretrain_batches(2L), a)
  expect_identical(batch_digests(a), ref_batches()$parallel)
  # each sample has its own stream, so the worker count doesn't matter
  expect_identical(capture_pretrain_batches(3L), a)
  for (x in a[lengths(lapply(a, dim)) == 4L]) {
    flat <- apply(x, 1L, paste, collapse = ",")
    expect_equal(anyDuplicated(flat), 0L)
  }
  # the caller's future plan is left alone
  expect_s3_class(future::plan(), "sequential")
})

test_that("workers get the static inputs once; batches carry only the specs", {
  skip_if_workers_differ()
  sizes <- new.env()
  sizes$call <- numeric(0); sizes$apply <- numeric(0)
  ser <- function(x) length(serialize(x, NULL))
  suppressMessages({
    trace("clusterCall", where = asNamespace("parallel"), print = FALSE,
          tracer = bquote(assign("call", c(get("call", envir = .(sizes)),
                                           .(ser)(list(fun, list(...)))),
                                 envir = .(sizes))))
    trace("clusterApply", where = asNamespace("parallel"), print = FALSE,
          tracer = bquote(assign("apply", c(get("apply", envir = .(sizes)),
                                            vapply(x, function(xi)
                                              .(ser)(list(fun, xi, list(...))),
                                              numeric(1))),
                                 envir = .(sizes))))
  })
  on.exit(suppressMessages({
    untrace("clusterCall", where = asNamespace("parallel"))
    untrace("clusterApply", where = asNamespace("parallel"))
  }), add = TRUE)
  # ~20 MB of val_real lives in pretrain_denoising_online()'s frame; under
  # the old dispatch every future carried it.
  val_real <- replicate(600, matrix(0.5, 128, 32), simplify = FALSE)
  capture_pretrain_batches(2L, val_real = val_real)
  cfg_bytes <- ser(make_toy_catalog())
  # one clusterCall carries the inputs (to every worker), once per run
  expect_equal(sum(sizes$call > cfg_bytes), 1L)
  expect_lt(max(sizes$call), 1e6)
  # per batch, each worker gets only its own specs + seeds: at most a few
  # H x W stripe / baseline matrices (32 KB each here), never the frame
  expect_true(length(sizes$apply) >= 6L)
  expect_lt(max(sizes$apply), 1e6)
  expect_lt(max(sizes$apply), ser(val_real) / 20)
})

test_that("workers are stopped when the run errors", {
  skip_on_cran()
  got <- new.env()
  suppressMessages(trace(".pretrain_start_workers",
                         where = asNamespace("classydms"), print = FALSE,
                         exit = bquote(assign("w", returnValue(),
                                              envir = .(got)))))
  on.exit(suppressMessages(untrace(".pretrain_start_workers",
                                   where = asNamespace("classydms"))),
          add = TRUE)
  # wrong val_real dims: fails after the workers are up
  expect_error(suppressWarnings(capture_pretrain_batches(
    2L, val_real = list(matrix(0, 3, 3)))))
  expect_false(is.null(got$w))
  expect_error(parallel::clusterCall(got$w$cl, Sys.getpid))
})

test_that("light contamination regime fills config defaults before scaling", {
  skip_on_cran()
  # Without n_per_sample / intensity_range, NULL * 0.3 used to give
  # numeric(0), so the generator's defaults never applied: NA intensities
  # and a "non-finite pixels" error.
  batches <- capture_pretrain_batches(
    1L, contamination_config = list(),
    contamination_schedule = list(light = 1))
  expect_true(length(batches) > 0L)
  for (x in batches) expect_true(all(is.finite(x)))
})
