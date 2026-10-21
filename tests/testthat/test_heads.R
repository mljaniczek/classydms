skip_if_not_installed("torch")
skip_if_not(torch::torch_is_installed())

toy_catalog <- function() {
  structure(list(
    compounds = tibble::tibble(compound_id = 1:5,
                               rt_loc = c(.1, .3, .5, .7, .9),
                               cv_loc = c(.2, .4, .5, .6, .8)),
    parameters = list(coord_mode = "fraction",
                      data_dims = cbind(rows = 1000, cols = 90))),
    class = "peak_catalog")
}

test_that("rank_auc matches known values", {
  expect_equal(rank_auc(c(0, 0, 1, 1), c(1, 2, 3, 4)), 1)
  expect_equal(rank_auc(c(0, 1), c(1, 1)), 0.5)
  expect_true(is.na(rank_auc(c(1, 1), c(1, 2))))
})

test_that("compound_head_geometry folds shift into prior width", {
  g0 <- compound_head_geometry(toy_catalog(), H = 1000, W = 100)
  g1 <- compound_head_geometry(toy_catalog(), H = 1000, W = 100, extra_rt_px = 20)
  expect_equal(g0$compound_id, 1:5)
  expect_true(all(g1$sd_rt_frac > g0$sd_rt_frac))
  expect_error(compound_head_geometry(toy_catalog(), 100, 10, compound_ids = 99))
})

test_that("presence head attends at catalog locations and round-trips", {
  geo <- compound_head_geometry(toy_catalog(), H = 1400, W = 90)
  enc <- dms_encoder()
  f <- enc(torch::torch_randn(2, 1, 1400, 90))
  ph <- dms_presence_head(64L, geo)
  expect_equal(ph(f)$shape, c(2, 5))
  attn <- ph$pool(f)$attn
  h <- f$size(3); w <- f$size(4)
  for (k in c(1, 5)) {
    am <- matrix(torch::as_array(attn[1, k, ]), nrow = h, byrow = TRUE)
    pk <- which(am == max(am), arr.ind = TRUE)[1, ]
    expect_lte(abs(pk[1] - geo$rt_frac[k] * h), 1.5)
    expect_lte(abs(pk[2] - geo$cv_frac[k] * w), 1.5)
  }
  tf <- tempfile(fileext = ".pt")
  torch::torch_save(ph, tf)
  expect_equal(torch::torch_load(tf)$compound_ids, 1:5)
})

test_that("classifier readouts produce logits", {
  geo <- compound_head_geometry(toy_catalog(), H = 256, W = 32)
  enc <- dms_encoder()
  ph <- dms_presence_head(64L, geo)
  x <- torch::torch_randn(3, 1, 256, 32)
  for (ro in c("gap", "attention", "presence")) {
    m <- build_pretrained_classifier(enc, readout = ro, presence_head = ph)
    m$eval()
    expect_equal(m(x)$shape, c(3, 2))
  }
  expect_error(build_pretrained_classifier(enc, readout = "presence"),
               "requires presence_head")
})
