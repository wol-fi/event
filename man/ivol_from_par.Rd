\name{ivol_from_par}
\alias{ivol_from_par}
\title{Implied volatility from model parameters via CF-based option pricing}
\description{
Computes Black--Scholes implied volatility by first pricing options under the model (via \code{opt_price})
and then inverting the Black--Scholes formula using a robust root-finding procedure.
}
\usage{
ivol_from_par(
  K,
  T = 0.9,
  par_x,
  par_y,
  r = 0,
  div = 0,
  type = c("call", "put"),
  S0 = NULL,
  vol_lower = 1e-10,
  vol_upper = 5,
  tol = 1e-10,
  maxiter = 200,
  ...
)
}
\arguments{
  \item{K}{Numeric vector of strike prices. Must be strictly positive.}
  \item{T}{Time to maturity in years. Must be strictly positive.}
  \item{par_x}{List of parameters for the \code{x}-block (passed to \code{opt_price} and possibly \code{get_S0}).}
  \item{par_y}{List of parameters for the \code{y}-block (passed to \code{opt_price} and possibly \code{get_S0}).}
  \item{r}{Continuously compounded risk-free rate.}
  \item{div}{Continuously compounded dividend yield (or foreign rate in FX setting).}
  \item{type}{Option type, either \code{"call"} or \code{"put"}.}
  \item{S0}{Spot price. If \code{NULL}, the function attempts to call \code{get_S0(par_x, par_y)}.}
  \item{vol_lower}{Lower bound for implied volatility search interval.}
  \item{vol_upper}{Initial upper bound for implied volatility search interval; may be expanded internally if needed.}
  \item{tol}{Absolute tolerance passed to \code{\link[stats]{uniroot}}.}
  \item{maxiter}{Maximum number of iterations passed to \code{\link[stats]{uniroot}}.}
  \item{...}{Additional arguments forwarded to \code{opt_price}.}
}
\details{
For each strike \code{K[i]}, the function:
\enumerate{
  \item Computes the model option price using \code{opt_price}.
  \item Checks no-arbitrage bounds for the corresponding Black--Scholes price.
  \item Solves for the Black--Scholes implied volatility using \code{\link[stats]{uniroot}}.
}
If the model price violates no-arbitrage bounds or the implied volatility cannot be bracketed,
\code{NA} is returned for that strike.
}
\value{
A numeric vector of implied volatilities of the same length as \code{K}. Entries may be \code{NA}.
}
\seealso{
\code{\link[stats]{uniroot}}
}
\examples{
\dontrun{
K <- 80:120
iv <- ivol_from_par(K, T = 30/365, par_x = par_x, par_y = par_y, r = 0, div = 0, type = "call")
plot(K, iv, type = "l")
}
}
\keyword{finance}
