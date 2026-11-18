# CDF derivatives retain their original evaluation points, order, and checks.
multiplierCdfJacobian <- function(points, F, theta, method.args, callback) {
  points <- as.matrix(points)
  Ftheta <- function(theta) F(points, theta)
  value <- if (is.null(callback)) numericalJacobian(Ftheta, theta, method.args)
           else callback(theta, points)
  if (!is.numeric(value) || !is.matrix(value) ||
      !identical(dim(value), c(NROW(points), length(theta))) ||
      any(!is.finite(value)))
    stop("CDF derivative must be a finite evaluation-point by parameter matrix")
  value
}

# One sorted cumulative-sum table serves both observation and grid statistics.
multiplierEcdf <- function(x, z) {
  ordering <- order(x)
  cumulative <- matrix(0, nrow(z), length(x) + 1L)
  for (i in seq_len(nrow(z))) cumulative[i, -1L] <- cumsum(z[i, ordering])
  list(x = x[ordering], cumulative = cumulative)
}

multiplierStat <- function(points, x, z, empirical, projection, F, theta,
                           method.args, callback, block.size = 128L) {
  points <- as.matrix(points)
  n <- NROW(x)
  derivative <- multiplierCdfJacobian(points, F, theta, method.args, callback)
  if (!is.null(empirical)) {
    # findInterval includes every tied observation <= the evaluation point.
    index <- findInterval(points[, 1], empirical$x) + 1L
    paths <- (empirical$cumulative[, index, drop = FALSE] -
      projection %*% t(derivative)) / n
    return(apply(paths, 1, function(w) c(max(abs(w)), n * mean(w * w))))
  }
  # Multivariate lower-orthant counts are computed in bounded column blocks.
  # Neither the full indicator matrix nor the full correction is constructed.
  ng <- NROW(points)
  maximum <- squares <- numeric(nrow(projection))
  for (first in seq.int(1L, ng, by = block.size)) {
    index <- seq.int(first, min(ng, first + block.size - 1L))
    bins <- vapply(index, function(j)
      rowSums(sweep(x, 2, points[j, ], "<=")) == NCOL(x), logical(n))
    paths <- (z %*% bins - projection %*% t(derivative[index, , drop = FALSE])) / n
    maximum <- pmax(maximum, apply(abs(paths), 1, max))
    squares <- squares + rowSums(paths * paths)
  }
  rbind(maximum, n * squares / ng)
}
