# Generic numerical differentiation for the default multiplier-bootstrap path.
# There is no distribution- or parameter-name-specific scaling.
derivative_failure <- function(status, parameter = NA_integer_, detail = "") {
  message <- paste0("Derivative selection: ", status,
    if (!is.na(parameter)) paste0(" for parameter ", parameter) else "", ". ", detail,
    " Supply explicit scales/steps or analytical derivatives; consider rescaling or reparameterizing.")
  stop(structure(list(message = message, call = NULL, status = status,
    parameter = parameter), class = c("fgof_derivative_failure", "error", "condition")))
}

stable_scale <- function(logdens, theta, target = 0.05, maxiter = 100L) {
  safe <- function(parm) tryCatch({
    value <- suppressWarnings(as.numeric(logdens(parm)))
    if (length(value) == 0L || any(!is.finite(value))) NULL else value
  }, error = function(e) NULL)
  baseline <- safe(theta)
  if (is.null(baseline)) derivative_failure("nonfinite_base")
  repeat.base <- safe(theta)
  if (!identical(baseline, repeat.base)) derivative_failure("nondeterministic_function")
  noise.floor <- 64 * .Machine$double.eps * max(1, sqrt(mean(baseline^2)))
  scales <- radii <- numeric(length(theta)); diagnostics <- list()
  for (j in seq_along(theta)) {
    h <- 1e-4 * if (theta[j] == 0) 1 else abs(theta[j])
    radius <- 0; last <- NULL; invalid.upper <- Inf; status <- "iteration_limit"
    for (iteration in seq_len(maxiter)) {
      plus <- minus <- theta
      plus[j] <- theta[j] + h; minus[j] <- theta[j] - h
      if (plus[j] == theta[j] || minus[j] == theta[j])
        derivative_failure("unrepresentable_pilot", j)
      fp <- safe(plus); fm <- safe(minus)
      valid <- !is.null(fp) && !is.null(fm) &&
        length(fp) == length(baseline) && length(fm) == length(baseline)
      if (!valid) {
        invalid.upper <- min(invalid.upper, h)
        if (!is.null(last) && last$signal < target / 2) {
          h <- last$h; signal <- last$signal; status <- "domain_limited"; break
        }
        h <- h / 2; next
      }
      radius <- max(radius, h)
      signal <- max(sqrt(mean((fp - baseline)^2)), sqrt(mean((fm - baseline)^2)))
      if (!is.finite(signal) || signal > 2 * target) h <- h / 2
      else if (signal < target / 2) {
        last <- list(h = h, signal = signal)
        if (2 * h >= invalid.upper) { status <- "domain_limited"; break }
        h <- h * 2
      } else { status <- "calibrated"; break }
    }
    if (status == "iteration_limit") derivative_failure(
      if (radius == 0) "no_two_sided_neighborhood" else
      if (!is.null(last) && last$signal <= noise.floor) "unresolved_pilot_signal" else status, j)
    if (!is.finite(signal) || signal <= noise.floor)
      derivative_failure("unresolved_pilot_signal", j)
    scales[j] <- h / signal
    radii[j] <- radius
    diagnostics[[j]] <- data.frame(parameter = j, pilot_step = h,
      rms_logdensity_change = signal, iterations = iteration, scale = scales[j],
      radius = radius, status = status)
  }
  list(scale = scales, radius = radii, diagnostics = do.call(rbind, diagnostics))
}

stable_jacobian <- function(func, theta, scale, cap = 0.1, method.args,
                            radius = rep(Inf, length(theta)),
                            rtol = 1e-5, atol = 1e-7, maxiter = 12L) {
  p <- length(theta)
  if (length(scale) != p || any(!is.finite(scale) | scale <= 0) ||
      length(radius) != p || any(is.na(radius) | radius <= 0))
    stop("scales/radii must be positive and match the parameter dimension")
  if (length(cap) != 1L || !is.finite(cap) || cap <= 0) stop("invalid relative cap")
  r <- method.args$r
  if (length(r) != 1L || r < 2 || r != as.integer(r) || method.args$v != 2)
    stop("Richardson differentiation requires integer r >= 2 and v = 2")
  safe <- function(parm) tryCatch({
    value <- suppressWarnings(as.numeric(func(parm)))
    if (!length(value) || any(!is.finite(value))) NULL else value
  }, error = function(e) NULL)
  baseline <- safe(theta)
  if (is.null(baseline)) derivative_failure("nonfinite_base")
  if (!identical(baseline, safe(theta))) derivative_failure("nondeterministic_function")
  initial <- abs(method.args$d * theta) +
    method.args$eps * (abs(theta) < method.args$zero.tol)
  ceiling <- pmin(cap * scale, radius)
  # Keep the historical answer only after an independent bounded check passes.
  original <- NULL
  if (length(initial) == p && all(initial > 0 & initial <= ceiling))
    original <- tryCatch(suppressWarnings(numDeriv::jacobian(func, theta,
      method.args = method.args)), error = function(e) NULL)
  value <- matrix(NA_real_, length(baseline), p)
  steps <- numeric(p); diagnostics <- list()
  for (j in seq_len(p)) {
    cache <- new.env(parent = emptyenv())
    assign(sprintf("%a", theta[j]), baseline, envir = cache)
    evaluate <- function(parm) {
      key <- sprintf("%a", parm)
      if (exists(key, cache, inherits = FALSE)) return(get(key, cache))
      point <- theta; point[j] <- parm
      out <- safe(point)
      if (!is.null(out) && length(out) != length(baseline)) out <- NULL
      assign(key, out, envir = cache)
      out
    }
    derivative <- function(h) {
      displacements <- h / 2^(0:(r - 1L))
      if (any(theta[j] + displacements == theta[j] |
              theta[j] - displacements == theta[j])) return(list(status = "resolution"))
      # Require the actual represented offsets to match the nominal stencil.
      mismatch <- max(abs((theta[j] + displacements - theta[j]) / displacements - 1),
        abs((theta[j] - (theta[j] - displacements)) / displacements - 1))
      if (mismatch > max(8 * .Machine$double.eps, rtol / 100)) return(list(status = "resolution"))
      one <- function(z) {
        out <- evaluate(theta[j] + h * z)
        if (is.null(out)) stop("invalid perturbation")
        out
      }
      args <- method.args; args$d <- 0; args$eps <- 1; args$zero.tol <- Inf
      out <- tryCatch(numDeriv::jacobian(one, 0, method.args = args), error = function(e) NULL)
      if (is.null(out) || any(!is.finite(out))) return(list(status = "domain"))
      list(status = "ok", value = as.numeric(out) / h)
    }
    discrepancy <- function(a, b) max(abs((a - b) * scale[j]) /
      (atol + rtol * pmax(abs(a * scale[j]), abs(b * scale[j]))))
    # Dyadic physical steps preserve the represented Richardson stencil at
    # large offsets; each candidate remains below both limits.
    h <- 2^floor(log2(ceiling[j]))
    if (!is.finite(h) || h <= 0) derivative_failure("unrepresentable_cap", j)
    previous <- older <- NULL; accepted <- FALSE; last.status <- "unstable"
    best <- Inf
    for (iteration in seq_len(maxiter)) {
      current <- derivative(h)
      last.status <- current$status
      if (current$status == "resolution") break
      if (current$status != "ok") { previous <- older <- NULL; h <- h / 2; next }
      if (!is.null(previous) && !is.null(older)) {
        err <- max(discrepancy(current$value, previous$value), discrepancy(previous$value, older$value))
        best <- min(best, err)
        if (err <= 1) { accepted <- TRUE; break }
      }
      older <- previous; previous <- current; h <- h / 2
    }
    if (!accepted) derivative_failure(if (last.status == "resolution")
      "insufficient_resolution_within_cap" else "unstable_within_cap", j)
    # The middle of three agreeing estimates avoids preferring the smallest
    # step merely because it was evaluated last.
    candidate <- previous$value
    selected.h <- h * 2
    retained <- !is.null(original) && all(is.finite(original[, j])) &&
      discrepancy(original[, j], candidate) <= 1
    value[, j] <- if (retained) original[, j] else candidate
    steps[j] <- if (retained) initial[j] else selected.h
    diagnostics[[j]] <- data.frame(parameter = j, step = steps[j], cap = ceiling[j],
      retained, iterations = iteration, discrepancy_ratio = best,
      status = if (retained) "validated_original" else "stable_dyadic")
  }
  list(value = value, steps = steps, changed = any(!vapply(diagnostics,
    function(z) z$retained, logical(1))), diagnostics = do.call(rbind, diagnostics))
}

# Supplied parameter scales still require a feasible density neighborhood.
# Their magnitude is user-controlled, so no sensitivity target is imposed.
derivativeRadius <- function(logdens, theta, scale, cap = 0.1, maxiter = 100L) {
  safe <- function(parm) tryCatch({
    value <- suppressWarnings(as.numeric(logdens(parm)))
    if (!length(value) || any(!is.finite(value))) NULL else value
  }, error = function(e) NULL)
  baseline <- safe(theta)
  if (is.null(baseline)) derivative_failure("nonfinite_base")
  if (!identical(baseline, safe(theta))) derivative_failure("nondeterministic_function")
  radius <- numeric(length(theta))
  for (j in seq_along(theta)) {
    h <- cap * scale[j]
    accepted <- FALSE
    for (iteration in seq_len(maxiter)) {
      plus <- minus <- theta
      plus[j] <- theta[j] + h; minus[j] <- theta[j] - h
      if (!is.finite(h) || h <= 0 || plus[j] == theta[j] || minus[j] == theta[j])
        derivative_failure("unrepresentable_pilot", j)
      fp <- safe(plus); fm <- safe(minus)
      if (!is.null(fp) && !is.null(fm) && length(fp) == length(baseline) &&
          length(fm) == length(baseline)) { accepted <- TRUE; break }
      h <- h / 2
    }
    if (!accepted) derivative_failure("no_two_sided_neighborhood", j)
    radius[j] <- h
  }
  radius
}

# One lazy context is shared by the score and each CDF evaluation set.
derivativeEvaluator <- function(logdens, parscale = NULL) {
  pilot <- current.theta <- NULL
  function(func, theta, method.args) {
    if (is.null(pilot) || !identical(theta, current.theta)) {
      if (is.null(parscale)) pilot <<- stable_scale(logdens, theta)
      else {
        if (!is.numeric(parscale) || length(parscale) != length(theta) ||
            any(!is.finite(parscale) | parscale == 0))
          stop("optimizer parscale must give a finite nonzero scale per parameter")
        scale <- abs(unname(parscale))
        pilot <<- list(scale = scale, radius = derivativeRadius(logdens, theta, scale))
      }
      current.theta <<- theta
    }
    stable_jacobian(func, theta, pilot$scale, cap = 0.1,
      method.args = method.args, radius = pilot$radius)$value
  }
}

numericalJacobian <- function(func, theta, method.args) {
  evaluator <- attr(method.args, "fgof.evaluator", exact = TRUE)
  if (is.null(evaluator)) numDeriv::jacobian(func, theta, method.args = method.args)
  else evaluator(func, theta, method.args)
}
