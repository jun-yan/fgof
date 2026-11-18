test_that("automatic derivatives work for custom models at difficult scales", {
  functions <- list(
    dfgofderiv = function(x, a, b, log = FALSE) dnorm(x, a, b, log = log),
    pfgofderiv = function(q, a, b) pnorm(q, a, b),
    rfgofderiv = function(n, a, b) rnorm(n, a, b))
  stopifnot(!any(vapply(names(functions), exists, logical(1),
    envir = .GlobalEnv, inherits = FALSE)))
  for (name in names(functions)) assign(name, functions[[name]], .GlobalEnv)
  on.exit(rm(list = names(functions), envir = .GlobalEnv))
  score <- function(theta, points) {
    z <- (points[, 1] - theta[1]) / theta[2]
    cbind(z / theta[2], (z^2 - 1) / theta[2])
  }
  cdf <- function(theta, points) {
    z <- (points[, 1] - theta[1]) / theta[2]
    cbind(-dnorm(z) / theta[2], -z * dnorm(z) / theta[2])
  }
  for (theta in list(c(0, 1e-6), c(1e6, 0.01), c(1e12, 1))) {
    set.seed(940009)
    x <- rnorm(80, theta[1], theta[2])
    args <- list(x = x, distr = "fgofderiv", nsim = 19, ng = 17,
      start = list(a = theta[1], b = theta[2]),
      optim.args = list(
        lower = c(theta[1] - 10 * theta[2], theta[2] / 100),
        upper = c(theta[1] + 10 * theta[2], theta[2] * 100),
        control = list(parscale = rep(theta[2], 2), maxit = 2000, factr = 1e4)))
    set.seed(950009)
    reference <- do.call(gof.test, c(args,
      list(derLogPdfWrtPar = score, derCdfWrtPar = cdf)))
    seed <- .Random.seed
    for (callbacks in list(list(), list(derLogPdfWrtPar = score),
                           list(derCdfWrtPar = cdf))) {
      set.seed(950009)
      numerical <- do.call(gof.test, c(args, callbacks))
      expect_identical(numerical$estimate, reference$estimate)
      expect_identical(numerical$statistic, reference$statistic)
      expect_lt(max(abs(numerical$stat.sim - reference$stat.sim)), 1e-6)
      expect_identical(numerical$p.value, reference$p.value)
      expect_identical(.Random.seed, seed)
    }
    # An explicit request retains the legacy numerical path, including failures.
    if (theta[1] == 0) expect_error(suppressWarnings(do.call(gof.test,
      c(args, list(numDeriv.args = fgof:::gradControl())))), "finite n by p")
  }
})

test_that("ordinary numerical results and parametric bootstrap remain compatible", {
  set.seed(941)
  x <- rnorm(60, 0.5, 1.2)
  for (simulation in c("mult", "pb")) {
    set.seed(942)
    automatic <- gof.test(x, "norm", 7, list(mean = 0.5, sd = 1.2),
      simulation = simulation, ng = 17)
    seed <- .Random.seed
    set.seed(942)
    explicit <- gof.test(x, "norm", 7, list(mean = 0.5, sd = 1.2),
      simulation = simulation, ng = 17, numDeriv.args = fgof:::gradControl())
    expect_identical(automatic, explicit)
    expect_identical(.Random.seed, seed)
  }
  h <- c(1e-4, 1e-8)
  controls <- list(d = 0, eps = h, zero.tol = Inf, r = 6)
  f <- function(theta) c(sin(theta[1]), theta[2]^2)
  expect_identical(fgof:::numericalJacobian(f, c(0, 1e-6), controls),
    numDeriv::jacobian(f, c(0, 1e-6), method.args = controls))
})

test_that("fixed spread and valid zero locations have automatic derivative scales", {
  set.seed(943)
  x <- rnorm(80, 0, 1e-6)
  args <- list(x = x, distr = "norm", nsim = 19, start = list(mean = 0),
    distr.args = list(sd = 1e-6), ng = 17,
    optim.args = list(lower = -1e-5, upper = 1e-5,
      control = list(parscale = 1e-6)))
  set.seed(944)
  reference <- do.call(gof.test, c(args, list(
    derLogPdfWrtPar = function(theta, points) matrix((points[, 1] - theta[1]) / 1e-12, ncol = 1),
    derCdfWrtPar = function(theta, points) matrix(-dnorm((points[, 1] - theta[1]) / 1e-6) / 1e-6, ncol = 1))))
  set.seed(944)
  automatic <- do.call(gof.test, args)
  expect_lt(max(abs(automatic$stat.sim - reference$stat.sim)), 1e-7)
  expect_identical(automatic$p.value, reference$p.value)
  # The production pilot also works with no optimizer scale supplied.
  f <- function(p) dnorm(x, p[1], 1e-6, log = TRUE)
  evaluator <- fgof:::derivativeEvaluator(f)
  value <- evaluator(function(p) pnorm(x, p[1], 1e-6), 0, fgof:::gradControl())
  expect_lt(max(abs(value[, 1] + dnorm(x / 1e-6) / 1e-6)) * 1e-6, 1e-7)
})

test_that("unstable, noisy, and boundary derivatives give actionable diagnostics", {
  x <- seq(-2, 2, length.out = 101)
  weak <- function(p) {
    if (p[1] <= 0) stop("outside domain")
    dnorm(x, p[1]^2, 1, log = TRUE)
  }
  pilot <- fgof:::stable_scale(weak, 1e-6)
  expect_identical(pilot$diagnostics$status, "domain_limited")
  expect_lt(pilot$radius, 1e-6)
  expect_error(fgof:::stable_jacobian(weak, 1e-6, pilot$scale,
    method.args = fgof:::gradControl(), radius = pilot$radius),
    "unstable_within_cap", class = "fgof_derivative_failure")
  boundary <- function(p) {
    if (p[1] < 0) stop("outside domain")
    dnorm(x, p[1], log = TRUE)
  }
  expect_error(fgof:::stable_scale(boundary, 0), "no_two_sided_neighborhood",
    class = "fgof_derivative_failure")
  set.seed(945)
  noisy <- function(p) pnorm(x, p[1]) + rnorm(length(x), sd = 1e-8)
  expect_error(fgof:::stable_jacobian(noisy, 0, 1,
    method.args = fgof:::gradControl()), "nondeterministic_function",
    class = "fgof_derivative_failure")
  expect_error(fgof:::stable_jacobian(function(p) pnorm(x, p[1]), 1e12, 1,
    cap = 1e-6, method.args = fgof:::gradControl()),
    "insufficient_resolution_within_cap", class = "fgof_derivative_failure")
})

test_that("correlated vector-parameter models retain analytical agreement", {
  sigma <- function(rho) matrix(c(1, rho, rho, 1), 2)
  functions <- list(
    dfgofcorr = function(x, param, log = FALSE)
      mvtnorm::dmvnorm(x, sigma = sigma(param[1]), log = log),
    pfgofcorr = function(q, param) apply(q, 1, function(point)
      as.numeric(mvtnorm::pmvnorm(upper = point, sigma = sigma(param[1]),
        algorithm = mvtnorm::TVPACK()))),
    rfgofcorr = function(n, param) mvtnorm::rmvnorm(n, sigma = sigma(param[1])))
  stopifnot(!any(vapply(names(functions), exists, logical(1),
    envir = .GlobalEnv, inherits = FALSE)))
  for (name in names(functions)) assign(name, functions[[name]], .GlobalEnv)
  on.exit(rm(list = names(functions), envir = .GlobalEnv))
  score <- function(theta, points) {
    rho <- theta[1]; cross <- points[, 1] * points[, 2]
    matrix(rho / (1 - rho^2) + ((1 + rho^2) * cross -
      rho * rowSums(points^2)) / (1 - rho^2)^2, ncol = 1)
  }
  cdf <- function(theta, points) matrix(mvtnorm::dmvnorm(points,
    sigma = sigma(theta[1])), ncol = 1)
  set.seed(946)
  x <- mvtnorm::rmvnorm(60, sigma = sigma(0.3))
  args <- list(x = x, distr = "fgofcorr", nsim = 9, start = list(param = 0.3),
    optim.args = list(lower = -0.8, upper = 0.8), gridStat = FALSE)
  set.seed(947)
  reference <- do.call(gof.test, c(args,
    list(derLogPdfWrtPar = score, derCdfWrtPar = cdf)))
  set.seed(947)
  automatic <- do.call(gof.test, args)
  expect_identical(automatic$estimate, reference$estimate)
  expect_lt(max(abs(automatic$stat.sim - reference$stat.sim)), 1e-6)
  expect_identical(automatic$p.value, reference$p.value)
})
