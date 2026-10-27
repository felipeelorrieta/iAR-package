library(iAR)

# ── iAR validators ────────────────────────────────────────────────────────────

test_that("iAR rejects non-scalar coef", {
  expect_error(iAR(times = 1:5, coef = c(0.9, 0)), "scalar")
})

test_that("iAR rejects coef <= 0 and suggests CiAR", {
  expect_error(iAR(times = 1:5, coef = -0.5), "CiAR")
})

test_that("iAR rejects coef >= 1", {
  expect_error(iAR(times = 1:5, coef = 1.0), "stationarity")
})

# ── CiAR / BiAR validators ────────────────────────────────────────────────────

test_that("CiAR rejects scalar coef", {
  expect_error(CiAR(times = 1:5, coef = 0.9), "length 2")
})

test_that("CiAR rejects coef with Mod >= 1", {
  expect_error(CiAR(times = 1:5, coef = c(0.8, 0.7)), "modulus")
})

test_that("BiAR rejects rho outside [-1, 1]", {
  expect_error(BiAR(times = 1:5, rho = 1.5), "\\[-1, 1\\]")
})

test_that("BiAR rejects coef with Mod >= 1", {
  expect_error(BiAR(times = 1:5, coef = c(0.9, 0.5)), "modulus")
})

# ── gentime ───────────────────────────────────────────────────────────────────

test_that("gentime rejects unknown distribution", {
  expect_error(gentime(NULL, n = 10, distribution = "poisson"), "distribution")
})

test_that("gentime NULL dispatch returns utilities object", {
  set.seed(1)
  o <- gentime(NULL, n = 10, distribution = "uniform")
  expect_s3_class(o, "iAR::utilities")
  expect_length(o@times, 10)
  expect_true(all(diff(o@times) > 0))
})

# ── kalman / loglik edge cases ─────────────────────────────────────────────────

test_that("kalman stops on series shorter than 3", {
  m <- iAR(family = "norm", times = 1:2, series = c(1, -1),
           series_esd = c(0, 0))
  expect_error(kalman(m), "at least 3")
})

test_that("kalman stops on constant series", {
  m <- iAR(family = "norm", times = 1:5, series = rep(1, 5),
           series_esd = rep(0, 5))
  expect_error(kalman(m), "constant")
})

# ── iAR-Gamma sim ─────────────────────────────────────────────────────────────

test_that("gamma sim stops when mean <= 0", {
  set.seed(1)
  o <- gentime(NULL, n = 20)
  m <- iAR(family = "gamma", times = o@times, coef = 0.8)
  m@mean     <- 0
  m@variance <- 1
  expect_error(sim(m), "strictly positive")
})

test_that("gamma sim works with non-integer mean and variance", {
  set.seed(42)
  o <- gentime(NULL, n = 50)
  m <- iAR(family = "gamma", times = o@times, coef = 0.8)
  m@mean     <- 1.5
  m@variance <- 0.5
  m <- sim(m)
  expect_length(m@series, 50)
  expect_true(all(m@series > 0))
})

# ── multidata dimension check ─────────────────────────────────────────────────

test_that("multidata rejects mismatched series_esd dimensions", {
  expect_error(
    multidata(times = 1:3,
              series     = matrix(1:6, nrow = 3, ncol = 2),
              series_esd = matrix(0,   nrow = 2, ncol = 3)),
    "Dimensions"
  )
})
