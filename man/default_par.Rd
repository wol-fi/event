\name{default_model_params}
\alias{default_X}
\alias{default_x}
\alias{default_x_bates}
\alias{default_x_heston}
\alias{default_x_BS}
\alias{default_Y}
\title{Default parameter lists for market-risk and event-risk blocks}
\description{
Constructs default parameter lists for the market-risk (X-block) and event-risk (Y-block) components.
User-supplied elements override defaults.
}
\usage{
default_X(par = list())

default_x(par = list())

default_x_bates(par = list())

default_x_heston(par = list())

default_x_BS(par = list())

default_Y(par = list())
}
\arguments{
  \item{par}{A named list of parameters to override defaults.}
}
\details{
\code{default_X()}, \code{default_x()} and \code{default_x_bates()} return Bates (Heston + lognormal jumps) defaults:
\itemize{
  \item \code{S0x = 100}
  \item \code{kappa = 1.5}, \code{theta = 0.04}, \code{sigma = 0.4}, \code{rho = -0.7}, \code{v0 = 0.04}
  \item \code{lambda = 0.2}, \code{muJ = -0.05}, \code{sigJ = 0.15}
}
\code{default_x_heston()} returns the same Heston parameters, but with jumps switched off
(\code{lambda = 0}, \code{muJ = 0}, \code{sigJ = 0}).

\code{default_x_BS()} returns Black--Scholes defaults in the same parameter format:
\code{kappa = 0}, \code{theta = 0}, \code{rho = 0}, \code{v0 = 0}, \code{lambda = 0}, \code{muJ = 0}, \code{sigJ = 0},
with \code{S0x = 100} and \code{sigma = 0.2} (interpreted as the Black--Scholes spot volatility).

\code{default_Y()} returns a list with event-risk parameters for the Wright--Fisher block, in particular
\code{p0} (initial/mean level) and \code{sigma_p} (volatility rate). In \code{cf_y} the Wright--Fisher
transition density is approximated by a Beta distribution parameterized by \code{p0} and \code{sigma_p}.
}
\value{
A named list of parameters. The returned list contains at least:
\itemize{
  \item \code{default_X/default_x/default_x_bates}: \code{S0x}, \code{kappa}, \code{theta}, \code{sigma}, \code{rho}, \code{v0}, \code{lambda}, \code{muJ}, \code{sigJ}
  \item \code{default_x_heston}: \code{S0x}, \code{kappa}, \code{theta}, \code{sigma}, \code{rho}, \code{v0}, \code{lambda}, \code{muJ}, \code{sigJ} (with jump parameters set to 0)
  \item \code{default_x_BS}: \code{S0x}, \code{kappa}, \code{theta}, \code{sigma}, \code{rho}, \code{v0}, \code{lambda}, \code{muJ}, \code{sigJ} (with non-BS parameters set to 0)
  \item \code{default_Y}: \code{p0}, \code{sigma_p}, \code{gamma}, \code{mh}, \code{sh}, \code{ml}, \code{sl}, \code{n_quad}
}
}
\examples{
par_x <- default_x()
par_x2 <- default_x(list(S0x=120, rho=-0.5))

par_h <- default_x_heston()
par_h2 <- default_x_heston(list(v0=0.03))

par_bs <- default_x_BS()
par_bs2 <- default_x_BS(list(sigma=0.25))

par_y <- default_Y()
par_y2 <- default_Y(list(p0=0.3))
}
