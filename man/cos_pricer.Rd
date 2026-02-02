\name{cos_pricer}
\alias{cos_pricer}
\title{COS option pricer for market-risk (X) times event-risk (Y)}
\description{
Prices European options using the COS method, where the underlying is modeled as
\eqn{S_T = X_T Y_T} with independent blocks \code{X} (market risk, Bates) and \code{Y} (event risk).
The characteristic function of \eqn{\log(S_T)} is built as \code{cf_x(u)*cf_y(u)}.
}
\usage{
cos_pricer(K, T, r, div = 0, par_x = list(), par_y = list(),
           type = c("call", "put"), N = 1024, L = 10)
}
\arguments{
  \item{K}{Numeric vector of strikes.}
  \item{T}{Time to maturity.}
  \item{r}{Continuously compounded risk-free rate.}
  \item{div}{Continuously compounded dividend yield (or foreign rate), default 0.}
  \item{par_x}{List of market-risk (Bates) parameters. Must include \code{S0x}.
See \code{\link{cf_x}} and \code{\link{default_X}}.}
  \item{par_y}{List of event-risk parameters.
See \code{\link{cf_y}} and \code{\link{default_Y}}.}
  \item{type}{Option type, \code{"call"} or \code{"put"}.}
  \item{N}{Number of COS terms (positive integer).}
  \item{L}{Truncation width parameter (positive).}
}
\details{
Let \eqn{Z = \log(S_T/K)}. The COS method approximates the discounted payoff expectation on a truncated interval
\eqn{[a,b]} determined from cumulants of \eqn{Z} (using \code{L}), and computes Fourier-cosine coefficients of the payoff.
The characteristic function used is
\deqn{\phi_{\log S}(u) = \phi_{\log X}(u)\,\phi_{\log Y}(u).}
For multiple strikes, the routine prices each strike using the same model inputs.
}
\value{
Numeric vector of option prices with the same length and order as \code{K}.
}
\seealso{
\code{\link{cf_x}}, \code{\link{cf_y}}, \code{\link{default_X}}, \code{\link{default_Y}}.
}
\examples{
par_x <- default_X(list(S0x=100))
par_y <- default_Y()

K <- 80:120
C <- cos_pricer(K, T=1, r=0, div=0, par_x=par_x, par_y=par_y, type="call")
}
