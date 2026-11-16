library(fgof)

set.seed(123)
x <- rlnorm(40, 1, 0.5)
start <- list(meanlog = 1, sdlog = 0.5)
set.seed(42)
default <- gof.test(x, "lnorm", 3, start, ng = 17)
set.seed(42)
explicit <- gof.test(x, "lnorm", 3, start, simulation = "mult",
                     qgrid = "empirical", ng = 17)
stopifnot(identical(default, explicit))

for (method in c("pb", "mult")) {
  set.seed(43)
  a <- gof.test(x, "lnorm", 3, start, simulation = method, ng = 17)
  set.seed(43)
  b <- gof.test(x, "lnorm", 3, start, simulation = method, ng = 17)
  stopifnot(identical(a, b), identical(dim(a$stat.sim), c(4L, 3L)),
            identical(names(a), c("statistic", "p.value", "estimate", "stat.sim")))
}

stopifnot(inherits(try(gof.test(x, "lnorm", 3, start, simulation = "typo"),
                       silent = TRUE), "try-error"),
          inherits(try(gof.test(x, "lnorm", 3, start, qgrid = "typo"),
                       silent = TRUE), "try-error"))

# Namespace calls must work when dependency packages are not attached.
detach("package:fgof")
detach("package:numDeriv")
detach("package:mvtnorm")
set.seed(42)
stopifnot(identical(default, fgof::gof.test(x, "lnorm", 3, start, ng = 17)))
