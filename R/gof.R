#### parse turns out to be quite time consuming
#### so better parse once and assign the expression to the body of a function
#### this is what fitdistr does in package MASS


#### get MLE
findMLE <- function(x, densfun, start, distr.args = list(), optim.args = list(), bootstrap = NULL) {
  if (missing(densfun) || !(is.function(densfun)) )
    stop("'densfun' must be supplied as a function")
  if (!is.list(distr.args) ||
      (length(distr.args) > 0L &&
       (is.null(names(distr.args)) || anyNA(names(distr.args)) ||
        any(names(distr.args) == "") || anyDuplicated(names(distr.args)))))
    stop("'distr.args' must be a list with unique, nonempty names")
  if (!is.list(optim.args) ||
      (length(optim.args) > 0L &&
       (is.null(names(optim.args)) || anyNA(names(optim.args)) ||
        any(names(optim.args) == "") || anyDuplicated(names(optim.args)))))
    stop("'optim.args' must be a list with unique, nonempty names")
  if (any(names(optim.args) %in% c("par", "fn")))
    stop("'optim.args' cannot set 'par' or 'fn'")
  if (any(names(distr.args) %in% names(start)))
    stop("'distr.args' cannot repeat parameters estimated in 'start'")
  log.arg <- "log" %in% names(formals(densfun))
  n <- NROW(x)

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
  dens <- function(parm) {
    parameters <- if (length(nm) == 1L) stats::setNames(list(parm), nm)
                  else stats::setNames(as.list(parm), nm)
    call.args <- c(list(x), parameters, distr.args)
    if (log.arg) call.args$log <- TRUE
    do.call(densfun, call.args)
  }

  objective <- if (log.arg) {
    function(parm) -sum(dens(parm))
  } else {
    function(parm) -sum(log(dens(parm)))
  }
  if (is.null(optim.args$method)) {
    optim.args$method <- if (any(c("lower", "upper") %in% names(optim.args)))
      "L-BFGS-B" else "Nelder-Mead"
  }
  run <- function(par, args) {
    tryCatch(do.call(stats::optim, c(list(par = par, fn = objective), args)),
             error = identity)
  }
  converged <- function(result) {
    !inherits(result, "error") && result$convergence == 0L &&
      all(is.finite(result$par)) && is.finite(result$value)
  }
  describe <- function(result) {
    if (inherits(result, "error")) return(conditionMessage(result))
    paste(c(paste("code", result$convergence), result$message), collapse = ", ")
  }
  res <- run(start, optim.args)
  detail <- NULL
  if (!is.null(bootstrap) && identical(optim.args$method, "L-BFGS-B") &&
      !inherits(res, "error") && res$convergence == 52L &&
      any(grepl("ABNORMAL_TERMINATION_IN_LNSRCH", res$message, fixed = TRUE)) &&
      is.null(optim.args$gr) && is.null(optim.args$control$ndeps) &&
      all(is.finite(res$par)) && is.finite(res$value)) {
    retry.args <- optim.args
    # ndeps is measured on the optimizer's parscale coordinate system.
    retry.args$control$ndeps <- rep(1e-5, length(start))
    retry <- run(res$par, retry.args)
    if (converged(retry)) {
      res <- retry
    } else {
      detail <- paste0("; retry: ", describe(retry))
    }
  }
  if (!converged(res)) {
    if (is.null(bootstrap) && inherits(res, "error")) stop(res)
    context <- if (!is.null(bootstrap))
      paste0(" in bootstrap replicate ", bootstrap) else ""
    stop("optimization failed", context, " (", optim.args$method, "; start = ",
         paste(format(unlist(start)), collapse = ", "), "): ", describe(res),
         detail, call. = FALSE)
  }
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
##' @param numDeriv.args Numerical Jacobian controls passed to \code{numDeriv}. Omission selects automatic step checks; explicit values are forwarded unchanged.
##' @param derCdfWrtPar Optional CDF derivative callback.
##' @param derLogPdfWrtPar Optional log-density derivative callback.
##' @param distr.args Named list of fixed arguments passed to the distribution functions.
##' @param optim.args Named list of arguments passed to \code{stats::optim}, such as \code{method}, \code{lower}, \code{upper}, and \code{control}.
##' @param refit.failure Action after a parametric bootstrap refit fails.
##'   Use \code{"stop"} (default) or \code{"report"} to attempt all replicates
##'   and return failure details. Ignored for the multiplier bootstrap.
##' @rawRd
##' \details{
##'   Fixed distribution parameters belong in \code{distr.args}; optimizer
##'   settings belong in \code{optim.args}. Both are named lists. For example,
##'   use \code{distr.args = list(sd = 2)} to fix a normal standard deviation
##'   and \code{optim.args = list(method = "BFGS")} to choose an optimizer
##'   method. Parameters listed in \code{start} are estimated and must not
##'   also appear in \code{distr.args}.
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
##'   \code{control$parscale} values supplied in \code{optim.args}. These are
##'   numerical scales,
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
##'   replicate. If L-BFGS-B reports an abnormal line-search termination,
##'   the refit is retried once on the same sample, starting from the
##'   optimizer's returned estimate when the estimate and objective are finite.
##'   The retry sets \code{ndeps = 1e-5} per parameter on the optimizer's
##'   \code{par / parscale} scale. Supplying \code{optim.args$control$ndeps}
##'   or an optimizer gradient via \code{optim.args$gr} skips the retry.
##'   A failed refit stops the
##'   test by default. With \code{refit.failure = "report"}, failed columns
##'   contain \code{NA}. If any refits fail, p-values are \code{NA} and the
##'   bounds treat every missing statistic first as a nonexceedance and then
##'   as an exceedance, using the full requested replicate count. The bounds
##'   describe this bootstrap run and are separate from Monte Carlo uncertainty.
##'   The multiplier bootstrap uses centered standard normal multipliers
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
##'   \item{refit.failures}{In parametric report mode, a list containing the
##'     failure count, fraction, replicate indices, and error messages.}
##'   \item{p.value.bounds}{In parametric report mode, a matrix of lower and
##'     upper bounds with one row per statistic. Both bounds equal the
##'     reported p-value when all refits succeed.}
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
                     numDeriv.args = gradControl(), ## for jacobian,
                     derCdfWrtPar = NULL,         ## function dF.dTheta(theta, x)
                     derLogPdfWrtPar = NULL,      ## function dlogf.dTheta(theta, x)
                     distr.args = list(),
                     optim.args = list(),
                     refit.failure = c("stop", "report")) {
#### x: the data vector to be tested
#### distr: character string for the name of the distribution
#### nsim: bootstrap sample size
#### start: a named list to be passed to fitdistr
  automatic.derivatives <- missing(numDeriv.args)
  simulation <- match.arg(simulation)
  report.failures <- simulation == "pb" && match.arg(refit.failure) == "report"
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
  if (!is.list(numDeriv.args))
    stop("'numDeriv.args' must be a list")
  if (!is.list(distr.args) ||
      (length(distr.args) > 0L &&
       (is.null(names(distr.args)) || anyNA(names(distr.args)) ||
        any(names(distr.args) == "") || anyDuplicated(names(distr.args)))))
    stop("'distr.args' must be a list with unique, nonempty names")
  if (!is.list(optim.args) ||
      (length(optim.args) > 0L &&
       (is.null(names(optim.args)) || anyNA(names(optim.args)) ||
        any(names(optim.args) == "") || anyDuplicated(names(optim.args)))))
    stop("'optim.args' must be a list with unique, nonempty names")
  if (any(names(optim.args) %in% c("par", "fn")))
    stop("'optim.args' cannot set 'par' or 'fn'")
  if (any(names(distr.args) %in% names(start)))
    stop("'distr.args' cannot repeat parameters estimated in 'start'")
  supported.optim.args <- setdiff(names(formals(stats::optim)), "...")
  invalid.optim.args <- setdiff(names(optim.args), supported.optim.args)
  if (length(invalid.optim.args))
    stop("unsupported 'optim.args': ", paste(invalid.optim.args, collapse = ", "))
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
  fit <- findMLE(x, ddistr, start, distr.args, optim.args)
  
  thetahat <- fit$estimate
  p <- length(thetahat)
  log.arg <- "log" %in% names(formals(ddistr))
  
  refit <- function(data, theta, bootstrap) {
    initial <- if (dim == 1L) as.list(theta) else list(param = theta)
    findMLE(data, ddistr, initial, distr.args, optim.args, bootstrap)$estimate
  }

  ## Forward only arguments accepted by each distribution function.
  dist_call_args <- function(fun, value, parm) {
    parameters <- if (dim == 1) as.list(parm) else list(parm)
    fixed <- distr.args[names(distr.args) %in% names(formals(fun)) |
                        "..." %in% names(formals(fun))]
    c(list(value), parameters, fixed)
  }
  pdistr.my <- function(x, parm) {
    do.call(pdistr, dist_call_args(pdistr, x, parm))
  }
  ddistr.my <- function(x, parm) {
    args <- dist_call_args(ddistr, x, parm)
    if (log.arg) args$log <- TRUE
    ret <- do.call(ddistr, args)
    if (log.arg) ret else log(ret)
  }
  rdistr.my <- function(n, parm) {
    args <- dist_call_args(rdistr, n, parm)
    if (dim > 1 && "dim" %in% names(formals(rdistr))) args$dim <- dim
    do.call(rdistr, args)
  }
  
  ## generate a grid for numerigal integration if gridStat == TRUE
  xg <- NULL
  if (gridStat) {
    ng1 <- ceiling(ng^(1 / dim))
    pgrid <- (1:ng1) / (ng1 + 1)
    if (dim == 1 && qgrid == "fitted") {
      xglist <- do.call(qdistr, dist_call_args(qdistr, pgrid, thetahat))
      xg <- as.matrix(expand.grid(xglist))
    }
    else {
      xglist <- list()
      for (i in 1:NCOL(x)) xglist[[i]] <- stats::quantile(x[,i], prob = pgrid)
      xg <- as.matrix(expand.grid(xglist))
    }
  }
  
  ## fxg <- ddistr.my(xg, thetahat) ## not needed
  make_empirical_cdf <- if (dim > 1) mecdf else stats::ecdf
  Fn <- make_empirical_cdf(x)
  Ftheta <- function (x) pdistr.my(x, thetahat)
  
  ## test statistics
  stat <- getStat(x, n, Fn, Ftheta)
  names(stat) <- c("ksa", "cvma")

  if (gridStat) {
    stat <- c(getStat(xg, n, Fn, Ftheta), stat)
    names(stat) <- c("ks", "cvm", "ksa", "cvma")
  }

  failures <- list()
  pb <- function(bootstrap) {
    x.sim <- rdistr.my(n, thetahat)
    theta.sim <- if (report.failures) {
      tryCatch(refit(x.sim, thetahat, bootstrap), error = identity)
    } else refit(x.sim, thetahat, bootstrap)
    if (inherits(theta.sim, "error")) {
      failures[[length(failures) + 1L]] <<- list(
        replicate = bootstrap, message = conditionMessage(theta.sim))
      return(stats::setNames(rep(NA_real_, length(stat)), names(stat)))
    }
    Fn.sim <- make_empirical_cdf(x.sim)
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
      (is.null(derLogPdfWrtPar) || is.null(derCdfWrtPar))) {
    parscale <- if (is.list(optim.args$control)) optim.args$control$parscale else NULL
    attr(numDeriv.args, "fgof.evaluator") <- derivativeEvaluator(logdens,
      parscale = parscale)
  }

  if (simulation == "pb") stat.sim <- vapply(seq_len(nsim), pb, stat)
  else {
    if (is.null(derLogPdfWrtPar))
      grad <- numericalJacobian(logdens, thetahat, numDeriv.args)
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
    stat.sim <- mult.sim(xg, x, pdistr.my, ddistr.my, thetahat, psi, nsim, gridStat = gridStat, method.args = numDeriv.args, derCdfWrtPar = derCdfWrtPar)
  }

  p.value <- (rowSums(stat.sim >= stat) + 0.5) / (nsim + 1)
  result <- list(statistic = stat, p.value = p.value, estimate = thetahat,
                 stat.sim = stat.sim)
  if (report.failures) {
    count <- length(failures)
    indices <- vapply(failures, function(failure) failure$replicate, integer(1))
    completed <- !seq_len(nsim) %in% indices
    exceeds <- rowSums(stat.sim[, completed, drop = FALSE] >= stat)
    result$refit.failures <- list(
      count = count, fraction = count / nsim,
      replicate = indices,
      message = vapply(failures, function(failure) failure$message, character(1)))
    result$p.value.bounds <- cbind(lower = (exceeds + 0.5) / (nsim + 1),
                                  upper = (exceeds + count + 0.5) / (nsim + 1))
  }
  result
}

mecdf <- function(x) {
  x <- as.matrix(x)
  function(points) copula::F.n(as.matrix(points), x, offset = 0, smoothing = "none")
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
