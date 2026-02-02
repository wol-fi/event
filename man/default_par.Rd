\name{default_model_params}
\alias{default_X}
\alias{default_Y}
\title{Default parameter lists for market-risk and event-risk blocks}
\description{
Constructs default parameter lists for the market-risk (X-block) and event-risk (Y-block) components.
User-supplied elements override defaults.
}
\usage{
default_X(par = list())

default_Y(par = list())
}
\arguments{
  \item{par}{A named list of parameters to override defaults.}
}
\details{
\code{default_X()} returns a list with Bates (Heston + jumps) parameters, including \code{S0x}.
\code{default_Y()} returns a list with event-risk parameters. If \code{p_var} is missing or not finite, it is set to
\code{0.6 * p0 * (1 - p0)}.
}
\value{
A named list of parameters. The returned list contains at least:
\itemize{
  \item \code{default_X}: \code{S0x}, \code{kappa}, \code{theta}, \code{sigma}, \code{rho}, \code{v0}, \code{lambda}, \code{muJ}, \code{sigJ}
  \item \code{default_Y}: \code{p0}, \code{p_var}, \code{gamma}, \code{mh}, \code{sh}, \code{ml}, \code{sl}, \code{n_quad}
}
}
\examples{
par_x <- default_X()
par_x2 <- default_X(list(S0x=120, rho=-0.5))

par_y <- default_Y()
par_y2 <- default_Y(list(p0=0.3))
}
