\name{cf_y}
\alias{cf_y}
\title{Event-risk characteristic function}
\description{
Computes the characteristic function of the log event-risk multiplier (Y-block) for a two-outcome event-risk model.
}
\usage{
cf_y(u, par_y)
}
\arguments{
  \item{u}{Numeric vector (or scalar). Points at which to evaluate the characteristic function.}
  \item{par_y}{List of event-risk parameters with elements:
    \describe{
      \item{p0}{Mean event probability in (0,1).}
      \item{p_var}{Variance of event probability. Must satisfy 0 < p_var < p0*(1-p0).}
      \item{gamma}{CRRA risk aversion parameter.}
      \item{mh}{Mean of log high-outcome multiplier.}
      \item{sh}{Std. dev. of log high-outcome multiplier.}
      \item{ml}{Mean of log low-outcome multiplier.}
      \item{sl}{Std. dev. of log low-outcome multiplier.}
      \item{n_quad}{Integer. Number of quadrature nodes (default: 64 if missing).}
    }
  }
}
\details{
The function uses a Gauss-Jacobi quadrature constructed via an eigen-decomposition of the tridiagonal Jacobi matrix.
Input validation:
\itemize{
  \item p0 must be in (0,1)
  \item p_var must satisfy 0 < p_var < p0*(1-p0)
}
}
\value{
Complex vector of the same length as \code{u}, giving \code{E[exp(1i*u*log(Y))]} under the pricing measure.
}
\examples{
par_y <- list(
  p0=0.4,
  p_var=0.6*0.4*(1-0.4),
  gamma=5,
  mh=log(1.10), sh=0.10,
  ml=log(0.90), sl=0.10,
  n_quad=64
)
u <- seq(0, 50, length.out=200)
phi <- cf_y(u, par_y)
}
