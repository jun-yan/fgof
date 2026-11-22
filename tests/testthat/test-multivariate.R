test_that("empirical lower-orthant counts handle ties and boundaries", {
  x <- rbind(c(0, 0), c(1, 0), c(1, 1), c(1, 1))
  points <- rbind(c(-1, -1), c(0, 0), c(1, 0), c(1, 1), c(2, 2))
  expect_identical(fgof:::mecdf(x)(points), c(0, 0.25, 0.5, 1, 1))
})

test_that("analytical t derivatives agree in dimensions two and three", {
  for (d in c(2L, 3L)) {
    set.seed(35)
    x <- matrix(rnorm(3 * d), ncol = d)
    theta <- c(rep(0, d), seq_len(d), rep(0.2, d * (d - 1) / 2))
    score <- numDeriv::jacobian(function(t)
      fgof:::dmymvt(x, t, "un", df = 5, log = TRUE), theta)
    cdf <- numDeriv::jacobian(function(t)
      fgof:::pmymvt(x, t, "un", df = 5), theta)
    expect_lt(max(abs(score - fgof:::gradLogPdfT(theta, x, df = 5))), 1e-6)
    expect_lt(max(abs(cdf - fgof:::gradCdfT(theta, x, df = 5))), 1e-6)
  }
})

test_that("both bivariate bootstrap paths return finite statistics", {
  set.seed(128)
  x <- matrix(rnorm(80), ncol = 2)
  for (method in c("mult", "pb")) {
    set.seed(34)
    out <- gof.test(x, "mymvnorm", 2, list(param = c(0, 0, 1, 1, 0)),
                    simulation = method,
                    distr.args = list(dispstr = "un"),
                    optim.args = list(lower = c(-5, -5, 0.1, 0.1, -0.9),
                                      upper = c(5, 5, 5, 5, 0.9)))
    empirical <- apply(x, 1, function(point) {
      mean(x[, 1] <= point[1] & x[, 2] <= point[2])
    })
    fitted <- fgof:::pmymvnorm(x, out$estimate, dispstr = "un")
    discrepancy <- empirical - fitted
    expect_equal(out$statistic,
                 c(ksa = max(abs(discrepancy)), cvma = sum(discrepancy^2)),
                 tolerance = 1e-6)
    expect_true(all(is.finite(unlist(out))))
    expect_identical(dim(out$stat.sim), c(2L, 2L))
  }
})
