test_that("fixed parameters reach statistics, fitted grids, and both bootstraps", {
  set.seed(124)
  x <- rnorm(60, 1, 2)
  for (method in c("pb", "mult")) {
    set.seed(31)
    a <- gof.test(x, "norm", 5, list(mean = 0), simulation = method,
                  gridStat = FALSE, sd = 2, method.args = list(), method = "BFGS")
    expected <- max(abs(ecdf(x)(x) - pnorm(x, a$estimate, sd = 2)))
    expect_lt(abs(a$estimate - mean(x)), 1e-5)
    expect_lt(abs(a$statistic["ksa"] - expected), 1e-12)
    expect_identical(names(a$p.value), names(a$statistic))
    expect_identical(rownames(a$stat.sim), names(a$statistic))

    y <- exp(x / 3)
    a <- gof.test(y, "lnorm", 3, list(sdlog = 1), simulation = method,
                  qgrid = "fitted", ng = 17, meanlog = 0, lower = 0.01)
    expect_lt(abs(a$estimate - sqrt(mean(log(y)^2))), 1e-4)
  }
})

test_that("gamma scale and rate names retain their parameterization", {
  set.seed(125)
  g <- rgamma(80, shape = 2, scale = 3)
  for (method in c("pb", "mult")) {
    set.seed(32)
    scale <- gof.test(g, "gamma", 3, list(shape = 2, scale = 3),
                      simulation = method, ng = 17, lower = c(0.01, 0.01))
    expect_lt(abs(mean(g) - prod(scale$estimate)), 1e-3)
    expected <- max(abs(ecdf(g)(g) - pgamma(g, shape = scale$estimate[1],
                                           scale = scale$estimate[2])))
    expect_lt(abs(scale$statistic["ksa"] - expected), 1e-12)
    set.seed(32)
    reversed <- gof.test(g, "gamma", 3, list(scale = 3, shape = 2),
                         simulation = method, ng = 17, lower = c(0.01, 0.01))
    expect_identical(scale, reversed)
    rate <- gof.test(g, "gamma", 3, list(shape = 2, rate = 1 / 3),
                     simulation = method, ng = 17, lower = c(0.01, 0.01))
    expect_lt(abs(rate$estimate[1] - scale$estimate[1]), 1e-3)
    expect_lt(abs(rate$estimate[2] - 1 / scale$estimate[2]), 1e-3)
    expect_lt(max(abs(rate$statistic - scale$statistic)), 1e-3)
  }
})

test_that("custom densities without a log argument retain likelihood behavior", {
  # Distribution lookup uses the search path, so install temporary test names.
  check_model <- function() {
    functions <- list(dfgofnolog = function(x, rate) dexp(x, rate),
                      pfgofnolog = function(q, rate) pexp(q, rate),
                      rfgofnolog = function(n, rate) rexp(n, rate))
    for (name in names(functions)) assign(name, functions[[name]], .GlobalEnv)
    on.exit(rm(list = names(functions), envir = .GlobalEnv))
    set.seed(129)
    z <- rexp(40, 2)
    for (method in c("pb", "mult")) {
      a <- gof.test(z, "fgofnolog", 3, list(rate = 2), simulation = method,
                    ng = 17, lower = 0.01)
      expect_lt(abs(a$estimate - 1 / mean(z)), 1e-4)
      expect_true(all(is.finite(unlist(a))))
    }
  }
  check_model()
})

test_that("analytical derivatives agree on distinct and equal-sized grids", {
  score <- function(theta, points) {
    z <- log(points[, 1]) - theta[1]
    cbind(z / theta[2]^2, -1 / theta[2] + z^2 / theta[2]^3)
  }
  cdf <- function(theta, points) {
    z <- (log(points[, 1]) - theta[1]) / theta[2]
    cbind(-dnorm(z) / theta[2], -z * dnorm(z) / theta[2])
  }
  set.seed(126)
  y <- rlnorm(40, 1, 0.5)
  for (grid in c("empirical", "fitted")) for (ng in c(17L, 40L)) {
    set.seed(33)
    numeric <- gof.test(y, "lnorm", 5, list(meanlog = 1, sdlog = 0.5),
                        ng = ng, qgrid = grid)
    set.seed(33)
    analytic <- gof.test(y, "lnorm", 5, list(meanlog = 1, sdlog = 0.5),
                         ng = ng, qgrid = grid, derCdfWrtPar = cdf,
                         derLogPdfWrtPar = score)
    expect_equal(numeric, analytic, tolerance = 1e-7)
  }
})
