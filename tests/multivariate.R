library(fgof)

# Independent lower-orthant counts, including ties and boundary points.
x <- rbind(c(0, 0), c(1, 0), c(1, 1), c(1, 1))
points <- rbind(c(-1, -1), c(0, 0), c(1, 0), c(1, 1), c(2, 2))
stopifnot(identical(fgof:::empiricalCdf(x)(points), c(0, 0.25, 0.5, 1, 1)))

# Preserve analytical t derivatives in both historically studied dimensions.
for (d in c(2L, 3L)) {
  set.seed(35)
  x <- matrix(rnorm(3 * d), ncol = d)
  theta <- c(rep(0, d), seq_len(d), rep(0.2, d * (d - 1) / 2))
  score <- numDeriv::jacobian(function(t)
    fgof:::dmymvt(x, t, "un", df = 5, log = TRUE), theta)
  cdf <- numDeriv::jacobian(function(t)
    fgof:::pmymvt(x, t, "un", df = 5), theta)
  stopifnot(max(abs(score - fgof:::gradLogPdfT(theta, x, df = 5))) < 1e-6,
            max(abs(cdf - fgof:::gradCdfT(theta, x, df = 5))) < 1e-6)
}

# End-to-end bivariate bootstrap paths with the internal normal wrappers.
set.seed(128)
x <- matrix(rnorm(80), ncol = 2)
for (method in c("mult", "pb")) {
  set.seed(34)
  out <- gof.test(x, "mymvnorm", 2, list(param = c(0, 0, 1, 1, 0)),
                  simulation = method, dispstr = "un",
                  lower = c(-5, -5, 0.1, 0.1, -0.9),
                  upper = c(5, 5, 5, 5, 0.9))
  stopifnot(all(is.finite(unlist(out))),
            identical(dim(out$stat.sim), c(2L, 2L)))
}
