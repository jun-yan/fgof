#### parse turns out to be quite time consuming
#### so better parse once and assign the expression to the body of a function
#### this is what fitdistr does in package MASS


#### get MLE
findMLE <- function(x, densfun, start, ...) {
  if (missing(densfun) || !(is.function(densfun)) )
    stop("'densfun' must be supplied as a function")
  log.arg <- "log" %in% names(formals(densfun))
  n <- NROW(x)
  
  myfn <- function(parm, ...) -sum(log(dens(parm, ...)))
  mylogfn <- function(parm, ...) -sum(dens(parm, ..., log = TRUE))

  ## correct the order of start to match the argument order in densfun
  nm <- names(start)
  f <- formals(densfun)
  args <- names(f)
  m <- match(nm, args)
  if (any(is.na(m))) 
    stop("'start' specifies names which are not arguments to 'densfun'")
  start <- start[order(m)]
  nm <- names(start)
  if (length(nm) == 1L) start <- c(unlist(start)) ## added for multivariate


  ## assuming for multivariate distribution, the parameters are stored
  ## in one parameter vector (e.g., param in dmymvnorm)
  dens <- function(parm, x, ...) {
    parameters <- if (length(nm) == 1L) stats::setNames(list(parm), nm)
                  else stats::setNames(as.list(parm), nm)
    do.call(densfun, c(list(x), parameters, list(...)))
  }
  
  Call <- match.call(expand.dots = TRUE)
  Call[[1L]] <- quote(stats::optim)
  Call$densfun <- Call$start <- Call$x <- NULL
  Call$fn <- if ("log" %in% args) mylogfn else myfn
  Call$par <- start
  Call$x <- x
  ## Call$control <- list(maxit=20000) ## should not be hard coded as this.

  ## Call$hessian <- TRUE
  if (is.null(Call$method)) {
    if (any(c("lower", "upper") %in% names(Call))) 
      Call$method <- "L-BFGS-B"
##    else if (length(start) > 1L) 
##      Call$method <- "BFGS"
##    else if (length(start) == 1) Call$method <- "Brent"
    else Call$method <- "Nelder-Mead" ## for dim > 1, this works better
  }

  res <- eval.parent(Call)
  if (res$convergence > 0L) 
    stop("optimization failed")
  ## sds <- sqrt(diag(solve(res$hessian)))
  list(estimate = res$par,
       ## sd = sds,
       loglik = -res$value, 
       n = n)
}

gradControl <- function(eps=1e-4, d=0.0001,
                        zero.tol=sqrt(.Machine$double.eps/7e-7),
                        r=6, v=2, show.details=FALSE) {
  list(eps=eps, d = d, zero.tol = zero.tol, r = r, v = v, show.details = show.details)
}

##' One Sample Goodness-of-fit Test
##'
##' Goodness-of-fit testing with multiplier or parametric bootstrap.
##' @param x A nonempty numeric vector or matrix of finite observations.
##' @param distr Name of the hypothesized distribution.
##' @param nsim Positive integer giving the number of bootstrap replicates.
##' @param start Nonempty named list of starting values for parameters to estimate.
##' @param simulation Bootstrap method: multiplier (`mult`) or parametric (`pb`).
##' @param ng Target number of grid points.
##' @param qgrid Grid construction method: empirical or fitted quantiles.
##' @param gridStat Whether to compute grid-based statistics in addition to observation-based statistics.
##' @param method.args Numerical Jacobian controls. Omission selects automatic step checks; explicit values are forwarded to numDeriv.
##' @param derCdfWrtPar Optional CDF derivative callback.
##' @param derLogPdfWrtPar Optional log-density derivative callback.
##' @param ... Fixed distribution arguments and optimizer arguments.
##' @rawRd
##' \details{
##'   Observations are assumed independent and identically distributed from a
##'   continuous distribution. Parameters in \code{start} are fitted by maximum
##'   likelihood. The optimizer defaults to Nelder-Mead, or L-BFGS-B when bounds
##'   are supplied. Choose appropriate starting values and bounds for the model;
##'   numerical derivatives may fail near parameter boundaries, and the multiplier
##'   procedure requires an invertible score covariance matrix.
##'
##'   Automatic numerical differentiation applies to standard and user-defined
##'   distributions alike. Parameter scales are inferred from two-sided changes
##'   in the fitted log density, or taken from the absolute optimizer
##'   \code{control$parscale} values when supplied. These are numerical scales,
##'   not parameter standard errors. A separate feasible radius records finite
##'   log-density probes; optimizer bounds do not by themselves define this radius.
##'   Candidate steps cannot exceed either that radius or 10 percent of the scale.
##'   Representable power-of-two steps are compared independently for log-density
##'   and CDF Jacobians. Three successive estimates must agree componentwise in
##'   scaled units, with absolute tolerance \code{1e-7} and relative tolerance
##'   \code{1e-5}. Richardson extrapolation retains six levels. A previous default
##'   derivative is kept when it meets these checks.
##'
##'   Selection stops with an actionable error if no reliable two-sided stencil
##'   is found or repeated function values reveal numerical randomness. It does
##'   not silently fall back to unchecked derivatives. Functions used for automatic
##'   selection must be deterministic and free of side effects. Approximate CDF
##'   integration, weak parameter sensitivity, constraints, and nonsmooth models
##'   can require analytical callbacks or explicitly selected numerical controls.
##'   Agreement between finite differences does not establish model regularity.
##'   No automatic selection is performed for parametric bootstrap. Analytical
##'   callbacks replace their corresponding derivatives; automatic checks continue
##'   for any remaining numerical Jacobian.
##'
##'   The parametric bootstrap simulates from the fitted model and refits every
##'   replicate. The multiplier bootstrap uses centered standard normal multipliers
##'   and an estimated parameter influence function. Grid-based statistics reuse
##'   the grid constructed from the original observations in both procedures.
##'
##'   The observation-based statistics \code{ksa} and \code{cvma} are respectively
##'   the maximum absolute empirical-minus-fitted CDF difference at observations
##'   and \eqn{n} times its mean square. The grid-based counterparts \code{ks} and
##'   \code{cvm} evaluate these expressions at empirical or fitted quantile points.
##'   The KS quantities use the unscaled differences, consistently in observed and
##'   simulated statistics. The empirical CDF in multiple dimensions counts
##'   componentwise inequalities.
##'
##'   Univariate models are the primary supported interface. Multivariate models
##'   require compatible distribution wrappers with a single vector parameter.
##'   The internal normal and t wrappers use integration algorithms for dimensions
##'   two and three; they do not provide general higher-dimensional support.
##'   Analytical t derivatives assume unrestricted correlations and fixed degrees
##'   of freedom. There is no interface for a test with all parameters known.
##'
##'   Set a random seed before calling the function for reproducible simulations
##'   within a fixed R and dependency environment. No seed is set internally.
##' }
##' @rawRd
##' \value{
##'   a list with the following components:
##'   \item{statistic}{Named vector of test statistics: \code{ks}, \code{cvm},
##'     \code{ksa}, \code{cvma} when \code{gridStat = TRUE}; otherwise \code{ksa},
##'     \code{cvma}.}
##'   \item{p.value}{Named vector of p-values in the same order as \code{statistic}.
##'     Computed as \code{(count + 0.5) / (nsim + 1)}, where \code{count} is the
##'     number of bootstrap statistics at least as large as the observed statistic.
##'     Values lie strictly between zero and one, including when no or all
##'     bootstrap statistics exceed the observed statistic.}
##'   \item{estimate}{vector of estimated parameter values}
##'   \item{stat.sim}{Matrix of bootstrap draws with one named row per statistic
##'     and \code{nsim} columns.}
##' }
##' @rawRd
##' \references{
##'   Kojadinovic I. and Yan J. (2012): Goodness-of-fit testing based on a weighted
##'   bootstrap: A fast large-sample alternative to the parametric
##'   bootstrap. Canadian Journal of Statistics, 40(3), 480--500.
##'   \doi{10.1002/cjs.11135}.
##' }
##' @rawRd
##' \author{
##'   Ivan Kojadinovic and Jun Yan.
##' }
##' @rawRd
##' \examples{
##' set.seed(123)
##' # Small simulation counts illustrate usage; use more for final inference.
##'
##' n <- 200
##' x <- rlnorm(200, 1, 1)
##'
##' system.time(m1 <- gof.test(x, "lnorm", 99, list(meanlog=1, sdlog=1), simulation="pb"))
##' system.time(m2 <- gof.test(x, "lnorm", 99, list(meanlog=1, sdlog=1), simulation="mult"))
##' apply(m1$stat.sim, 1, summary)
##' apply(m2$stat.sim, 1, summary)
##'
##' system.time(m1 <- gof.test(x, "gamma", 99, list(shape=1, rate=1), simulation="pb"))
##' system.time(m2 <- gof.test(x, "gamma", 99, list(shape=1, rate=1), simulation="mult"))
##' apply(m1$stat.sim, 1, summary)
##' apply(m2$stat.sim, 1, summary)
##'
##'
##' system.time(m1 <- gof.test(x, "exp", 99, list(rate=1), simulation="pb"))
##' system.time(m2 <- gof.test(x, "exp", 99, list(rate=1), simulation="mult"))
##' apply(m1$stat.sim, 1, summary)
##' apply(m2$stat.sim, 1, summary)
##'
##' }
##' @rawRd
##' \keyword{ htest }
##' @export
gof.test <- function(x, distr, nsim, start,
                     simulation = c("mult", "pb"),
                     ng = 1000,
                     qgrid = c("empirical", "fitted"), ## only matters for dim = 1
                     gridStat = NULL,  ## for CvM statistics, dim >= 2
                     method.args = gradControl(), ## for jacobian,
                     derCdfWrtPar = NULL,         ## function dF.dTheta(theta, x)
                     derLogPdfWrtPar = NULL,      ## function dlogf.dTheta(theta, x)
                     ...) {
#### x: the data vector to be tested
#### distr: character string for the name of the distribution
#### nsim: bootstrap sample size
#### start: a named list to be passed to fitdistr
  automatic.derivatives <- missing(method.args)
  simulation <- match.arg(simulation)
  qgrid <- match.arg(qgrid)
  if (!is.numeric(x) || !(is.null(dim(x)) || is.matrix(x)) ||
      NROW(x) < 1L || NCOL(x) < 1L || any(!is.finite(x)))
    stop("'x' must be a nonempty numeric vector or matrix of finite values")
  for (arg in c("nsim", "ng")) {
    value <- get(arg)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value) || value > .Machine$integer.max)
      stop("'", arg, "' must be a positive integer")
  }
  if (!is.list(start) || length(start) == 0L ||
      is.null(names(start)) || anyNA(names(start)) || any(names(start) == "") ||
      anyDuplicated(names(start)) ||
      any(!vapply(start, function(z) is.numeric(z) && length(z) > 0L &&
                   all(is.finite(z)), logical(1))))
    stop("'start' must be a nonempty named list of finite numeric parameters")
  if (NCOL(x) == 1L && any(lengths(start) != 1L))
    stop("each univariate parameter in 'start' must be a scalar")
  if (!is.null(gridStat) &&
      (length(gridStat) != 1L || is.na(gridStat) ||
       !(is.logical(gridStat) || is.numeric(gridStat)) ||
       !gridStat %in% c(FALSE, TRUE)))
    stop("'gridStat' must be TRUE or FALSE")
  for (fun in list(derCdfWrtPar, derLogPdfWrtPar))
    if (!is.null(fun) && !is.function(fun))
      stop("derivative callbacks must be functions or NULL")
  if (!is.list(method.args))
    stop("'method.args' must be a list; supply it explicitly when passing optim's 'method'")
  dots <- list(...)
  if (length(dots) > 0L && (is.null(names(dots)) || any(names(dots) == "")))
    stop("arguments in '...' must be named")
  n <- NROW(x)
  dim <- NCOL(x)
  if (is.null(gridStat)) gridStat <- if (dim == 1) TRUE else FALSE

  rdistr <- get(paste("r", distr, sep=""), mode = "function")
  pdistr <- get(paste("p", distr, sep=""), mode = "function")
  ddistr <- get(paste("d", distr, sep=""), mode = "function")
  if (dim == 1 && qgrid == "fitted") ## only for dim = 1
    qdistr <- get(paste("q", distr, sep=""), mode = "function") 
  
  if (!is.matrix(x)) x <- as.matrix(x)
  ## estimation
  fit <- findMLE(x, ddistr, start, ...)
  
  thetahat <- fit$estimate
  p <- length(thetahat)
  log.arg <- "log" %in% names(formals(ddistr))
  
  findMLE.call <- if (dim == 1) parse(text = paste("findMLE(x, ddistr, start = as.list(theta), ...)")) else parse(text = paste("findMLE(x, ddistr, start = list(param=theta), ...)"))

  ## Keep estimated parameters named, and forward fixed distribution arguments.
  fixed.args <- list(...)
  fixed.args <- fixed.args[!names(fixed.args) %in% names(formals(stats::optim))]
  distribution.args <- function(fun, value, parm) {
    parameters <- if (dim == 1) as.list(parm) else list(parm)
    fixed <- fixed.args[names(fixed.args) %in% names(formals(fun)) |
                        "..." %in% names(formals(fun))]
    c(list(value), parameters, fixed)
  }
  pdistr.my <- function(x, parm) {
    do.call(pdistr, distribution.args(pdistr, x, parm))
  }
  ddistr.my <- function(x, parm) {
    args <- distribution.args(ddistr, x, parm)
    if (log.arg) args$log <- TRUE
    ret <- do.call(ddistr, args)
    if (log.arg) ret else log(ret)
  }
  rdistr.my <- function(n, parm) {
    args <- distribution.args(rdistr, n, parm)
    if (dim > 1 && "dim" %in% names(formals(rdistr))) args$dim <- dim
    do.call(rdistr, args)
  }
  
  ## generate a grid for numerigal integration if gridStat == TRUE
  xg <- NULL
  if (gridStat) {
    ng1 <- ceiling(ng^(1 / dim))
    pgrid <- (1:ng1) / (ng1 + 1)
    if (dim == 1 && qgrid == "fitted") {
      xglist <- do.call(qdistr, distribution.args(qdistr, pgrid, thetahat))
      xg <- as.matrix(expand.grid(xglist))
    }
    else {
      xglist <- list()
      for (i in 1:NCOL(x)) xglist[[i]] <- stats::quantile(x[,i], prob = pgrid)
      xg <- as.matrix(expand.grid(xglist))
    }
  }
  
  ## fxg <- ddistr.my(xg, thetahat) ## not needed
  if (dim > 1) ecdf <- empiricalCdf
  Fn <- stats::ecdf(x)
  Ftheta <- function (x) pdistr.my(x, thetahat)
  
  ## test statistics
  stat <- getStat(x, n, Fn, Ftheta)
  names(stat) <- c("ksa", "cvma")

  if (gridStat) {
    stat <- c(getStat(xg, n, Fn, Ftheta), stat)
    names(stat) <- c("ks", "cvm", "ksa", "cvma")
  }

  pb <- function() {
    x.sim <- rdistr.my(n, thetahat)
    theta.sim <- eval(findMLE.call, list(theta = thetahat, x = x.sim))$estimate
    Fn.sim <- stats::ecdf(x.sim)
    Ftheta.sim <- function(x) pdistr.my(x, theta.sim)
    stat <- getStat(x.sim, n, Fn.sim, Ftheta.sim)
    if (gridStat) {
      stat <- c(getStat(xg, n, Fn.sim, Ftheta.sim), stat)
    }
    names(stat) <- if (gridStat) c("ks", "cvm", "ksa", "cvma")
                   else c("ksa", "cvma")
    stat
  }
  logdens <- function(theta) ddistr.my(x, theta)  
  llk <- function(theta) sum(logdens(theta))
  if (automatic.derivatives && simulation == "mult" &&
      (is.null(derLogPdfWrtPar) || is.null(derCdfWrtPar)))
    attr(method.args, "fgof.evaluator") <- derivativeEvaluator(logdens,
      parscale = dots$control$parscale)

  if (simulation == "pb") stat.sim <- replicate(nsim, pb())
  else {
    if (is.null(derLogPdfWrtPar))
      grad <- numericalJacobian(logdens, thetahat, method.args)
    else grad <- derLogPdfWrtPar(thetahat, x)
    ## a matrix of n by p
    ## hinv <- solve( - hessian(llk, theta) / n ) ## a matrix p by p
    ## if (fisher == "G") hinv <- solve(var(t(grad))) else hinv <- solve( - hessian(llk, theta) / n )
    ## hinv <- solve(var(grad))
    ## psi <- hinv %*% t(grad)
    if (!is.numeric(grad) || !is.matrix(grad) || !identical(dim(grad), c(n, p)) ||
        any(!is.finite(grad)))
      stop("log-density derivative must be a finite n by p matrix")
    psi <- tryCatch(solve(stats::var(grad), t(grad)), error = function(e)
      stop("cannot invert score covariance: ", conditionMessage(e), call. = FALSE))
    stat.sim <- mult.sim(xg, x, pdistr.my, ddistr.my, thetahat, psi, nsim, gridStat = gridStat, method.args = method.args, derCdfWrtPar = derCdfWrtPar)
  }

  p.value <- (rowSums(stat.sim >= stat) + 0.5) / (nsim + 1)
  list(statistic = stat, p.value = p.value, estimate = thetahat, stat.sim = stat.sim)
}

empiricalCdf <- function(x) {
  x <- as.matrix(x)
  function(points) {
    points <- as.matrix(points)
    vapply(seq_len(NROW(points)), function(i)
      mean(rowSums(sweep(x, 2, points[i, ], "<=")) == NCOL(x)), numeric(1))
  }
}

influFun <- function(xg, x, psi, F, theta, method.args, derCdfWrtPar) {
  ## xg: the grid point at which the function is to be evaluated
  ## x: observed data, n by d
  ## psi: influ of mle, a matrix of p by n
  ## F: cdf of X
  xg <- as.matrix(xg)
  x <- as.matrix(x)
  n <- NROW(x)
  d <- NCOL(x)
  ng <- NROW(xg)

  Ftheta <- function(theta) F(xg, theta)
  if (is.null(derCdfWrtPar))
    dF.dtheta <- numericalJacobian(Ftheta, theta, method.args)
  else dF.dtheta <- derCdfWrtPar(theta, xg)
  if (!is.numeric(dF.dtheta) || !is.matrix(dF.dtheta) ||
      !identical(dim(dF.dtheta), c(ng, length(theta))) ||
      any(!is.finite(dF.dtheta)))
    stop("CDF derivative must be a finite evaluation-point by parameter matrix")
  
  idx1 <- rep(1:n, ng)
  idx2 <- rep(1:ng, each = n)
  mat1 <- x[idx1,,drop = FALSE]
  mat2 <- xg[idx2,,drop = FALSE]
  bins <- matrix(ifelse(rowSums(mat1 <= mat2) == d, 1, 0), n, ng)
  ## influ.1 <- bins - matrix(Ftheta(theta), n, ng, byrow = TRUE) ## matrix of n by ng
  ## the second part matrix(...) is not needed when zcenter == TRUE
  influ.1 <- bins
  influ.2 <- t(dF.dtheta %*% psi)
  return(influ.1 - influ.2) ## a matrix of n by ng
}


influStat <- function(infl, z) {
  n <- ncol(z)
  simobj <- z %*% infl / n
  stat <- apply(simobj, 1, function(w) c(max(abs(w)), n * mean(w * w)))
  stat
}

mult.sim <- function(xg, x, pdistr.my, ddistr.my, thetahat, psi, nsim,
                     ## zcenter = TRUE, ## should always be TRUE
                     gridStat = TRUE,
                     method.args, derCdfWrtPar) {
  n <- NROW(x)
  d <- NCOL(x)
  z <- matrix(stats::rnorm(n * nsim), nrow = n)
  z <- t(scale(z, center = TRUE, scale = FALSE))

  # Associativity avoids a dense n-by-ng parameter-correction matrix.
  projection <- z %*% t(psi)
  empirical <- if (d == 1L) multiplierEcdf(x[, 1], z) else NULL
  if (d == 1L) z <- NULL
  stat.sim <- multiplierStat(x, x, z, empirical, projection, pdistr.my,
    thetahat, method.args, derCdfWrtPar)
  rownames(stat.sim) <- c("ksa", "cvma")
  if (gridStat) {
    stat.sim <- rbind(multiplierStat(xg, x, z, empirical, projection, pdistr.my,
      thetahat, method.args, derCdfWrtPar), stat.sim)
    rownames(stat.sim) <- c("ks", "cvm", "ksa", "cvma")
  }
  stat.sim
}

getStat <- function(x, n, Fn, Ftheta) {
  emp <- Fn(x) - Ftheta(x)
  ## n <- NROW(x)
  ##  c(ks = max(abs(emp)), cvm = n * weighted.mean(emp * emp, fxg),
  ##  wrong! because xg is not uniform grid
  c(ks = max(abs(emp)), cvm = n * mean(emp * emp))
}
