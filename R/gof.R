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
    parameters <- if (length(nm) == 1L) setNames(list(parm), nm)
                  else setNames(as.list(parm), nm)
    do.call(densfun, c(list(x), parameters, list(...)))
  }
  
  Call <- match.call(expand.dots = TRUE)
  Call[[1L]] <- as.name("optim")
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
  fixed.args <- fixed.args[!names(fixed.args) %in% names(formals(optim))]
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
      for (i in 1:NCOL(x)) xglist[[i]] <- quantile(x[,i], prob = pgrid)
      xg <- as.matrix(expand.grid(xglist))
    }
  }
  
  ## fxg <- ddistr.my(xg, thetahat) ## not needed
  if (dim > 1) ecdf <- empiricalCdf
  Fn <- ecdf(x)
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
    Fn.sim <- ecdf(x.sim)
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
    psi <- tryCatch(solve(var(grad), t(grad)), error = function(e)
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
  p <- length(thetahat)
  n <- NROW(x)
  d <- NCOL(x)
  ## xg is null unless gridStat == TRUE
  
  z <- matrix(rnorm(n * nsim), nrow = n)
  ## if (zcenter) z <- t(scale(z, center = TRUE, scale = FALSE)) ## nsim by n
  ## else z <- t(z)  ## nsim by n
  z <- t(scale(z, center = TRUE, scale = FALSE)) ## nsim by n

  infla <- influFun(x, x, psi, pdistr.my, thetahat, method.args, derCdfWrtPar)
  stat.sim <- influStat(infla, z)
  rownames(stat.sim) <- c("ksa", "cvma")
  if (gridStat) {
    infl <- influFun(xg, x, psi, pdistr.my, thetahat, method.args, derCdfWrtPar)
    stat.sim <- rbind(influStat(infl, z), stat.sim)
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
