test_that("invalid inputs, callback shapes, and optimization failures are clear", {
  set.seed(127)
  x <- rlnorm(30)
  args <- list(x = x, distr = "lnorm", nsim = 3,
               start = list(meanlog = 0, sdlog = 1), ng = 17)
  fails <- function(change, pattern) {
    a <- args
    a[names(change)] <- change
    expect_error(do.call(gof.test, a), pattern)
  }
  for (bad in list(numeric(), c(x, NA), c(x, Inf), letters[1:3]))
    fails(list(x = bad), "finite values")
  for (arg in c("nsim", "ng")) for (bad in list(0, -1, 1.5, NA, Inf, 1:2))
    fails(setNames(list(bad), arg), "positive integer")
  for (bad in list(list(), list(0, 1), list(meanlog = NA),
                   list(meanlog = c(0, 1)), list(meanlog = 0, meanlog = 1)))
    fails(list(start = bad), "start")
  for (bad in list(NA, c(TRUE, FALSE), "TRUE", 2))
    fails(list(gridStat = bad), "gridStat")
  fails(list(derCdfWrtPar = 1), "callbacks")
  fails(list(derLogPdfWrtPar = function(theta, x) matrix(0, 1, 1)),
        "log-density derivative")
  fails(list(derCdfWrtPar = function(theta, x) matrix(0, 1, 1)), "CDF derivative")
  fails(list(x = rep(1, 30), distr = "exp", start = list(rate = 1),
             numDeriv.args = list(),
             optim.args = list(method = "Brent", lower = 0.01, upper = 10)),
        "score covariance")
  fails(list(optim.args = list(control = list(maxit = 1))), "optimization failed")
  fails(list(distr.args = c(sd = 1, sd = 2)), "unique, nonempty names")
  fails(list(optim.args = list(unknown = 1)), "unsupported 'optim.args'")
  fails(list(optim.args = list(fn = identity)), "cannot set 'par' or 'fn'")
  fails(list(distr.args = list(meanlog = 0)), "cannot repeat parameters")
})

test_that("one replicate preserves output names and dimensions for both grids", {
  set.seed(127)
  x <- rlnorm(30)
  for (method in c("mult", "pb")) for (grid in c(FALSE, TRUE)) {
    out <- gof.test(x, "lnorm", 1, list(meanlog = 0, sdlog = 1),
                    simulation = method, gridStat = grid, ng = 17)
    expected <- if (grid) c("ks", "cvm", "ksa", "cvma") else c("ksa", "cvma")
    expect_identical(dim(out$stat.sim), c(length(expected), 1L))
    expect_true(all(out$p.value > 0 & out$p.value < 1))
    expect_identical(names(out$statistic), expected)
    expect_identical(names(out$p.value), expected)
    expect_identical(rownames(out$stat.sim), expected)
  }
})
