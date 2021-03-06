% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fwelnet.R
\name{fwelnet}
\alias{fwelnet}
\title{Fit feature-weighted elastic net}
\usage{
fwelnet(
  x,
  y,
  z,
  lambda = NULL,
  family = c("gaussian", "binomial"),
  alpha = 1,
  standardize = TRUE,
  max_iter = 1,
  ave_mode = 1,
  thresh_mode = 1,
  t = 1,
  a = 0.5,
  thresh = 1e-04,
  verbose = FALSE
)
}
\arguments{
\item{x}{Input matrix, of dimension \code{nobs x nvars}; each row is
an observation vector.}

\item{y}{Response variable. Quantitative for \code{family = "gaussian"}. For
\code{family="binomial"}, should be a numeric vector consisting of 0s and 1s.}

\item{z}{Feature of features matrix, with dimension \code{nvars x nfeaturevars}.}

\item{lambda}{A user supplied \code{lambda} sequence. Typical usage is to
have the program compute its own \code{lambda} sequence; supplying a value of
lambda overrides this.}

\item{family}{Response type. Either \code{"gaussian"} (default) for linear
regression or \code{"binomial"} for logistic regression.}

\item{alpha}{The elastic net mixing hyperparameter, a real value number
between 0 and 1 (inclusive). Default value is 1.}

\item{standardize}{If \code{TRUE}, the columns of the input matrix are
standardized before the algorithm is run. Default is \code{TRUE}.}

\item{max_iter}{The number of iterations for the optimization. Default is 1.}

\item{ave_mode}{If equal to 1 (default), the gradient descent direction for
\code{theta} is the mean gradient across the lambda values. If equal to 2,
it is the component-wise median gradient across the lambda values.}

\item{thresh_mode}{If equal to 1 (default), backtracking line search for
\code{theta} is done so that the mean objective function (across lambda
values) decreases. If equal to 2, it is done so that the median objective
function decreases.}

\item{t}{The initial step size for \code{theta} backtracking line search
(default value is 1).}

\item{a}{The factor by which step size is decreased in \code{theta}
backtracking line search (default value is 0.5).}

\item{thresh}{If the mean/median objective function does not decrease by at
least this factor, we terminate the optimization early. Default is 1e-4.}

\item{verbose}{If \code{TRUE}, prints information to console as model is
being fit. Default is \code{FALSE}.}
}
\value{
An object of class \code{"fwelnet"}.
\item{beta}{A \code{p x length(lambda)} matrix of coefficients.}
\item{theta}{Theta value, a \code{nfeaturevars x 1} matrix.}
\item{a0}{Intercept sequence of length \code{length(lambda)}.}
\item{lambda}{The actual sequence of \code{lambda} values used.}
\item{nzero}{The number of non-zero coefficients for each value of
\code{lambda}.}
\item{family}{Response type.}
\item{call}{The call that produced this object.}
\item{obj}{A \code{max_iter + 1} by \code{length(lambda)} matrix of
objective function values. (Number of rows could be fewer if the optimization
stopped early.)}
}
\description{
Fit a model with feature-weighted elastic net for a path of lambda values.
Fits linear and logistic regression models.
}
\details{
\code{fwelnet} always mean centers the columns of the \code{x} matrix. If
\code{standardize=TRUE}, \code{fwelnet} will also scale the columns to have
standard deviation 1. In all cases, the \code{beta} coefficients returned are
for the original \code{x} values (i.e. uncentered and unscaled).
}
\examples{
set.seed(1)
n <- 100; p <- 20
x <- matrix(rnorm(n * p), n, p)
beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
y <- x \%*\% beta + rnorm(n)
z <- cbind(1, abs(beta) + rnorm(p))

fwelnet(x, y, z)
fwelnet(x, y, z, ave_mode = 2)
fwelnet(x, y, z, ave_mode = 2, thresh_mode = 2)

}
