\name{cf_x}
\alias{cf_x}
\title{Market-risk (Bates) characteristic function}
\description{
Computes the characteristic function of \eqn{\log(X_T)} under a Bates model (Heston stochastic volatility with Merton-style lognormal jumps).
}
\usage{
cf_x(u, T, r, div = 0, par_x)
}
\arguments{
  \item{u}{Numeric vector (or scalar). Points at which to evaluate the characteristic function.}
  \item{T}{Time to maturity.}
  \item{r}{Continuously compounded risk-free rate.}
  \item{div}{Continuously compounded dividend yield (or foreign rate), default 0.}
  \item{par_x}{List of Bates parameters with elements:
    \describe{
      \item{S0x}{Initial level \eqn{X_0} (spot for the X-block).}
      \item{kappa}{Mean-reversion speed of variance.}
      \item{theta}{Long-run mean of variance.}
      \item{sigma}{Volatility of variance.}
      \item{rho}{Correlation between Brownian motions.}
      \item{v0}{Initial variance.}
      \item{lambda}{Jump intensity.}
      \item{muJ}{Mean of jump size in log-space.}
      \item{sigJ}{Std. dev. of jump size in log-space.}
    }
  }
}
\details{
Returns the Bates characteristic function
\deqn{\phi_X(u) = \mathbb{E}\left[\exp\{ i u \log(X_T)\}\right].}
The drift adjustment uses \eqn{\kappa_J = \exp(\mu_J + 0.5 \sigma_J^2) - 1}.
}
\value{
Complex vector of the same length as \code{u}, giving \code{E[exp(1i*u*log(X_T))]} under the pricing measure.
}
\examples{
par_x <- list(
  S0x=100,
  kappa=1.5, theta=0.04, sigma=0.4, rho=-0.7, v0=0.04,
  lambda=0.2, muJ=-0.05, sigJ=0.15
)
u <- seq(0, 50, length.out=200)
phi <- cf_x(u, T=1, r=0.01, div=0, par_x=par_x)
}
