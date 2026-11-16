test_that("default choices preserve explicit multiplier behavior", {
  set.seed(123)
  x <- rlnorm(40, 1, 0.5)
  start <- list(meanlog = 1, sdlog = 0.5)
  set.seed(42)
  default <- gof.test(x, "lnorm", 3, start, ng = 17)
  set.seed(42)
  explicit <- gof.test(x, "lnorm", 3, start, simulation = "mult",
                       qgrid = "empirical", ng = 17)
  expect_identical(default, explicit)
  expect_error(gof.test(x, "lnorm", 3, start, simulation = "typo"), "arg")
  expect_error(gof.test(x, "lnorm", 3, start, qgrid = "typo"), "arg")
})

test_that("both bootstrap methods reproduce seeded results and return layout", {
  set.seed(123)
  x <- rlnorm(40, 1, 0.5)
  start <- list(meanlog = 1, sdlog = 0.5)
  for (method in c("pb", "mult")) {
    set.seed(43)
    a <- gof.test(x, "lnorm", 3, start, simulation = method, ng = 17)
    set.seed(43)
    b <- gof.test(x, "lnorm", 3, start, simulation = method, ng = 17)
    expect_identical(a, b)
    expect_identical(dim(a$stat.sim), c(4L, 3L))
    expect_identical(names(a), c("statistic", "p.value", "estimate", "stat.sim"))
  }
  expect_identical(names(formals(gof.test)),
    c("x", "distr", "nsim", "start", "simulation", "ng", "qgrid", "gridStat",
      "method.args", "derCdfWrtPar", "derLogPdfWrtPar", "..."))
})
