\name{initial_levels}
\alias{get_Y0}
\alias{get_S0}
\title{Initial levels implied by the event-risk block}
\description{
Convenience functions to compute the initial event-risk multiplier \code{Y0} and the combined initial level \code{S0 = S0x * Y0}.
}
\usage{
get_Y0(par_y = list())

get_S0(par_y = list(), par_x = list())
}
\arguments{
  \item{par_y}{List of event-risk parameters passed to \code{cf_y} (or to the event-risk CF factory, if used).}
  \item{par_x}{List of market-risk parameters containing at least \code{S0x}.}
}
\details{
\code{get_Y0()} computes \eqn{Y_0 = \mathbb{E}[Y]} using the characteristic function of \eqn{\log(Y)}:
\deqn{Y_0 = \Re\{\phi_Y(-i)\}, \quad \phi_Y(u) = \mathbb{E}[\exp(i u \log Y)].}
\code{get_S0()} returns \code{S0x * get_Y0(par_y)}.

Both functions are intended as utilities for model parameterization and consistency checks.
}
\value{
\code{get_Y0()} returns a numeric scalar \code{Y0}.
\code{get_S0()} returns a numeric scalar \code{S0}.
}
\examples{
par_y <- default_Y()
par_x <- default_X(list(S0x=100))

Y0 <- get_Y0(par_y)
S0 <- get_S0(par_y, par_x)
}
