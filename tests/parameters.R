library(fgof)

set.seed(124)
x <- rnorm(60, 1, 2)
for (method in c("pb", "mult")) {
  set.seed(31)
  a <- gof.test(x, "norm", 5, list(mean = 0), simulation = method,
                gridStat = FALSE, sd = 2, method.args = list(), method = "BFGS")
  expected <- max(abs(ecdf(x)(x) - pnorm(x, a$estimate, sd = 2)))
  stopifnot(abs(a$estimate - mean(x)) < 1e-5,
            abs(a$statistic["ksa"] - expected) < 1e-12,
            identical(names(a$p.value), names(a$statistic)),
            identical(rownames(a$stat.sim), names(a$statistic)))

  # Check a fixed argument in a non-first parameter position and fitted grid.
  y <- exp(x / 3)
  a <- gof.test(y, "lnorm", 3, list(sdlog = 1), simulation = method,
                qgrid = "fitted", ng = 17, meanlog = 0, lower = 0.01)
  stopifnot(abs(a$estimate - sqrt(mean(log(y)^2))) < 1e-4)
}

set.seed(125)
g <- rgamma(80, shape = 2, scale = 3)
for (method in c("pb", "mult")) {
  set.seed(32)
  scale <- gof.test(g, "gamma", 3, list(shape = 2, scale = 3),
                    simulation = method, ng = 17, lower = c(0.01, 0.01))
  stopifnot(abs(mean(g) - prod(scale$estimate)) < 1e-3,
            abs(scale$statistic["ksa"] - max(abs(ecdf(g)(g) -
              pgamma(g, shape = scale$estimate[1], scale = scale$estimate[2])))) < 1e-12)
  set.seed(32)
  reversed <- gof.test(g, "gamma", 3, list(scale = 3, shape = 2),
                       simulation = method, ng = 17, lower = c(0.01, 0.01))
  stopifnot(identical(scale, reversed))
  rate <- gof.test(g, "gamma", 3, list(shape = 2, rate = 1 / 3),
                   simulation = method, ng = 17, lower = c(0.01, 0.01))
  stopifnot(abs(rate$estimate[1] - scale$estimate[1]) < 1e-3,
            abs(rate$estimate[2] - 1 / scale$estimate[2]) < 1e-3,
            max(abs(rate$statistic - scale$statistic)) < 1e-3)
}

# Custom densities without a log argument retain the non-log likelihood path.
dnolog <- function(x, rate) dexp(x, rate)
pnolog <- function(q, rate) pexp(q, rate)
rnolog <- function(n, rate) rexp(n, rate)
set.seed(129)
z <- rexp(40, 2)
for (method in c("pb", "mult")) {
  a <- gof.test(z, "nolog", 3, list(rate = 2), simulation = method,
                ng = 17, lower = 0.01)
  stopifnot(abs(a$estimate - 1 / mean(z)) < 1e-4,
            all(is.finite(unlist(a))))
}

# Numerical and analytical lognormal derivatives, including distinct grids.
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
  stopifnot(isTRUE(all.equal(numeric, analytic, tolerance = 1e-7)))
}
