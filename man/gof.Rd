\name{gof.test}
\alias{gof.test}
\title{
One Sample Goodness-of-fit Test
}
\description{
  Goodness-of-fit test for one sample data with two bootstrap methods:
  multiplier bootstrap and parametric bootstrap. The multipler bootstrap
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
    a vector or matrix of data to be tested.
  }
  \item{distr}{
    name of the hypothesized distribution, assuming that {dpr} functions
    are available for density, distribution, and random number
    generation. For example, if "foo" is the hypothesized distribution,
    then functions "dfoo", "pfoo", and "rfoo" must be available.
  }
  \item{nsim}{
    the bootstrap sample size.
  }
  \item{start}{
    a list of named components for the starting value of the parameter estimate.
  }
  \item{ng}{
    The size of the grid points for numerical integration.
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
    Not used for one dimensional distribution now.
  }
  \item{method.args}{
    a list of arguments to control the "method.args" for \code{jocobian}
    from package numDeriv. The defaul is:
    eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7),
    r=6, v=2, show.details=FALSE
  }
  \item{derCdfWrtPar}{
    derivative function of the cdf with respect to parameters.
    The dedault is NULL, in which case, numerical derivative is used.
  }
  \item{derLogPdfWrtPar}{
    derivative function of the logPdf with respect to parameters.
    The dedault is NULL, in which case, numerical derivative is used.
  }
  \item{\dots}{
   other arguments to be passed to {dpr} function of the distribution. 
  }
}
\details{
  The function has been fully tested with dimension 1. For higher
  dimensions, wrapper functions are often needed. Examples will be
  provided in future releases.
}
\value{
  a list with the following components:
  \item{statistic}{vector of test statistics}
  \item{p.value}{vector of p-values of the test statistics}
  \item{estimate}{vector of estimated parameter values}
  \item{stat.sim}{matrix of bootstrap draws of the test statistics}
}
\references{
  Kojadinovic I. and Yan J.: Goodness-of-fit testing based on a weighted
  bootstrap: A fast large-sample alternative to the parametric
  bootstrap. Canadian Journal of Statistics. Forthcoming.
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
