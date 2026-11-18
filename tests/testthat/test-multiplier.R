test_that("weighted empirical sums preserve ties, boundaries, and point order", {
  x <- c(1, 0, 1, -1, 2, 0)
  z <- rbind(c(0.5, -0.5, 1, -1, 0.25, -0.25),
             c(-1, 1, -0.5, 0.5, -0.25, 0.25))
  points <- c(Inf, 0, -Inf, 1, 2, -1, 0, 0.5)
  state <- fgof:::multiplierEcdf(x, z)
  actual <- state$cumulative[, findInterval(points, state$x) + 1L, drop = FALSE]
  expected <- vapply(points, function(point) rowSums(z[, x <= point, drop = FALSE]),
    numeric(nrow(z)))
  expect_identical(actual, expected)
})

test_that("factored and blocked statistics match the dense influence definition", {
  set.seed(68001)
  for (d in 1:3) for (nsim in c(1L, 7L)) {
    x <- matrix(round(rnorm(31 * d), 1), ncol = d)
    points <- matrix(round(rnorm(257 * d), 1), ncol = d)
    points[1, ] <- Inf; points[2, ] <- -Inf
    z <- t(scale(matrix(rnorm(NROW(x) * nsim), nrow = NROW(x)),
      center = TRUE, scale = FALSE))
    psi <- matrix(rnorm(2 * NROW(x)), 2)
    theta <- c(a = 0, b = 1)
    callback <- function(theta, points) cbind(plogis(points[, 1]), plogis(points[, d]))
    F <- function(points, theta) stop("analytical callback should bypass numerical CDF")
    dense <- fgof:::influFun(points, x, psi, F, theta, list(), callback)
    expected <- fgof:::influStat(dense, z)
    empirical <- if (d == 1L) fgof:::multiplierEcdf(x[, 1], z) else NULL
    actual <- fgof:::multiplierStat(points, x, z, empirical, z %*% t(psi),
      F, theta, list(), callback)
    expect_equal(unname(actual), unname(expected), tolerance = 1e-11)
    expect_identical(dim(actual), c(2L, nsim))
    one <- fgof:::multiplierStat(points[3, , drop = FALSE], x, z, empirical,
      z %*% t(psi), F, theta, list(), callback)
    expected.one <- fgof:::influStat(fgof:::influFun(points[3, , drop = FALSE],
      x, psi, F, theta, list(), callback), z)
    expect_equal(unname(one), unname(expected.one), tolerance = 1e-11)
  }
})
