\name{initial_levels}
\alias{get_Y0}
\alias{get_S0}
\title{Initial levels implied by the event-risk block}
\description{
Convenience functions to compute the initial event-risk multiplier \code{Y0} and the combined initial level \code{S0 = S0x * Y0}.
}
\usage{
get_Y0(par_y)

get_S0(par_x, par_y)
}
\arguments{
  \item{par_y}{List of event-risk parameters passed to \code{make_cf_y} (or \code{cf_y}, depending on your implementation). Must contain \code{tau}.}
  \item{par_x}{List of market-risk parameters containing at least \code{S0x}.}
}
\details{
\code{get_Y0()} computes \eqn{Y_0 = \mathbb{E}[Y]} using the characteristic function of \eqn{\log(Y)}:
\deqn{Y_0 = \Re\{\phi_{\log Y}(-i)\}, \quad \phi_{\log Y}(u) = \mathbb{E}[\exp(i u \log Y)].}
In code this is evaluated as \code{Re(make_cf_y(tau, par_y)(-1i))}.

\code{get_S0()} returns \code{par_x$S0x * get_Y0(par_y)}.

Both functions are intended as utilities for model parameterization and consistency checks.
}
\value{
\code{get_Y0()} returns a numeric scalar \code{Y0}.
\code{get_S0()} returns a numeric scalar \code{S0}.
}
\examples{
\dontrun{
par_y <- default_y()
par_x <- default_x(list(S0x = 100))

Y0 <- get_Y0(par_y)
S0 <- get_S0(par_x, par_y)
}
}
\keyword{models}
\keyword{math}
