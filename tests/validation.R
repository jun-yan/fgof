library(fgof)

set.seed(127)
x <- rlnorm(30)
args <- list(x = x, distr = "lnorm", nsim = 3,
             start = list(meanlog = 0, sdlog = 1), ng = 17)
fails <- function(change, pattern) {
  a <- args
  a[names(change)] <- change
  error <- tryCatch(do.call(gof.test, a), error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
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
           method.args = list(), method = "Brent", lower = 0.01, upper = 10), "score covariance")
fails(list(control = list(maxit = 1)), "optimization failed")

for (method in c("mult", "pb")) for (grid in c(FALSE, TRUE)) {
  out <- gof.test(x, "lnorm", 1, args$start, simulation = method,
                  gridStat = grid, ng = 17)
  expected <- if (grid) c("ks", "cvm", "ksa", "cvma") else c("ksa", "cvma")
  stopifnot(identical(dim(out$stat.sim), c(length(expected), 1L)),
            identical(names(out$statistic), expected),
            identical(names(out$p.value), expected),
            identical(rownames(out$stat.sim), expected))
}
