\name{gof.test}
\alias{gof.test}
\title{
One Sample Goodness-of-fit Test
}
\description{
  Goodness-of-fit test for one sample data with two bootstrap methods:
  multiplier bootstrap and parametric bootstrap. The multiplier bootstrap
  method is faster, using generic numeric derivatives.
}
\usage{
gof.test(x, distr, nsim, start, 
         simulation = c("mult", "pb"),
         ng = 1000, qgrid = c("empirical", "fitted"),
         gridStat = NULL, method.args = gradControl(),
         derCdfWrtPar = NULL, derLogPdfWrtPar = NULL, ...)
}
\arguments{
  \item{x}{
    a nonempty numeric vector or matrix of finite observations.
  }
  \item{distr}{
    name of the hypothesized distribution, assuming that \code{d}, \code{p}, and \code{r} functions
    are available for density, distribution, and random number
    generation. For example, if "foo" is the hypothesized distribution,
    then functions "dfoo", "pfoo", and "rfoo" must be available.
    A quantile function "qfoo" is also required for a fitted univariate grid.
  }
  \item{nsim}{
    the number of bootstrap replicates, a positive integer.
  }
  \item{start}{
    a nonempty list of named components for the starting value of the parameter
    estimate. Univariate components must be finite numeric scalars named as in
    the density function. All parameters in this list are estimated; it does
    not specify a fully known null distribution. Other distribution parameters
    can be held fixed through \code{...}.
  }
  \item{ng}{
    A positive integer specifying the target number of grid points. In dimension
    \eqn{d}, the Cartesian grid has \code{ceiling(ng^(1 / d))^d} points.
  }
  \item{simulation}{
    method of bootstrap: "mult" for multiplier bootstrap; "pb" for
    parametric bootstrap.
  }
  \item{qgrid}{
    method to construct the grid for numerical integration:
    "empirical" for empirical quantiles; "fitted" for fitted quantiles.
  }
  \item{gridStat}{
    Whether to compute grid-based statistics as well as observation-based
    statistics. The default \code{NULL} selects \code{TRUE} in dimension one
    and \code{FALSE} in higher dimensions.
  }
  \item{method.args}{
    a list of arguments to control \code{jacobian}
    from package numDeriv. The default is:
    eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7),
    r=6, v=2, show.details=FALSE
  }
  \item{derCdfWrtPar}{
    derivative function of the CDF with respect to parameters. It is called as
    \code{derCdfWrtPar(theta, points)} and must return a finite numeric matrix
    with one row per evaluation point and one column per estimated parameter.
    The default is \code{NULL}, in which case numerical derivatives are used.
  }
  \item{derLogPdfWrtPar}{
    derivative function of the log density with respect to parameters. It is
    called as \code{derLogPdfWrtPar(theta, x)} and must return a finite numeric
    matrix with one row per observation and one column per estimated parameter.
    Both callbacks receive a numeric parameter vector and a matrix of points,
    and must incorporate any fixed distribution parameters themselves.
    The default is \code{NULL}, in which case numerical derivatives are used.
  }
  \item{\dots}{
    Named fixed distribution arguments and arguments to \code{optim}, such as
    \code{lower}, \code{upper}, and \code{control}. Distribution arguments are
    forwarded to the corresponding density, CDF, quantile, and generator when
    accepted by their formal arguments or \code{...}. To pass \code{method}
    to \code{optim}, also supply \code{method.args} explicitly (for example,
    \code{method.args = list(), method = "BFGS"}) to prevent partial matching
    to \code{method.args}.
  }
}
\details{
  Observations are assumed independent and identically distributed from a
  continuous distribution. Parameters in \code{start} are fitted by maximum
  likelihood. The optimizer defaults to Nelder-Mead, or L-BFGS-B when bounds
  are supplied. Choose appropriate starting values and bounds for the model;
  numerical derivatives may fail near parameter boundaries, and the multiplier
  procedure requires an invertible score covariance matrix.

  The parametric bootstrap simulates from the fitted model and refits every
  replicate. The multiplier bootstrap uses centered standard normal multipliers
  and an estimated parameter influence function. Grid-based statistics reuse
  the grid constructed from the original observations in both procedures.

  The observation-based statistics \code{ksa} and \code{cvma} are respectively
  the maximum absolute empirical-minus-fitted CDF difference at observations
  and \eqn{n} times its mean square. The grid-based counterparts \code{ks} and
  \code{cvm} evaluate these expressions at empirical or fitted quantile points.
  The KS quantities use the unscaled differences, consistently in observed and
  simulated statistics. The empirical CDF in multiple dimensions counts
  componentwise inequalities.

  Univariate models are the primary supported interface. Multivariate models
  require compatible distribution wrappers with a single vector parameter.
  The internal normal and t wrappers use integration algorithms for dimensions
  two and three; they do not provide general higher-dimensional support.
  Analytical t derivatives assume unrestricted correlations and fixed degrees
  of freedom. There is no interface for a test with all parameters known.

  Set a random seed before calling the function for reproducible simulations
  within a fixed R and dependency environment. No seed is set internally.
}
\value{
  a list with the following components:
  \item{statistic}{Named vector of test statistics: \code{ks}, \code{cvm},
    \code{ksa}, \code{cvma} when \code{gridStat = TRUE}; otherwise \code{ksa},
    \code{cvma}.}
  \item{p.value}{Named vector of p-values in the same order as \code{statistic}.
    The archived convention divides the number of bootstrap statistics at least
    as large as the observed statistic by \code{nsim + 1}.}
  \item{estimate}{vector of estimated parameter values}
  \item{stat.sim}{Matrix of bootstrap draws with one named row per statistic
    and \code{nsim} columns.}
}
\references{
  Kojadinovic I. and Yan J. (2012): Goodness-of-fit testing based on a weighted
  bootstrap: A fast large-sample alternative to the parametric
  bootstrap. Canadian Journal of Statistics, 40(3), 480--500.
  \doi{10.1002/cjs.11135}.
}
\author{
  Ivan Kojadinovic and Jun Yan.
}
\examples{
set.seed(123)

n <- 200
x <- rlnorm(200, 1, 1)

system.time(m1 <- gof.test(x, "lnorm", 1000, list(meanlog=1, sdlog=1), simulation="pb"))
system.time(m2 <- gof.test(x, "lnorm", 200, list(meanlog=1, sdlog=1), simulation="mult"))
apply(m1$stat.sim, 1, summary)
apply(m2$stat.sim, 1, summary)

system.time(m1 <- gof.test(x, "gamma", 1000, list(shape=1, rate=1), simulation="pb"))
system.time(m2 <- gof.test(x, "gamma", 200, list(shape=1, rate=1), simulation="mult"))
apply(m1$stat.sim, 1, summary)
apply(m2$stat.sim, 1, summary)


system.time(m1 <- gof.test(x, "exp", 1000, list(rate=1), simulation="pb"))
system.time(m2 <- gof.test(x, "exp", 200, list(rate=1), simulation="mult"))
apply(m1$stat.sim, 1, summary)
apply(m2$stat.sim, 1, summary)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
