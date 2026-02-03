\name{cos_price}
\alias{cos_price}
\title{European option pricing via the COS method}
\description{
Prices European call or put options using the Fourier COS method applied to the
characteristic function of the log-price. The log-price CF is formed as the product
of two components via \code{make_cf_x()} and \code{make_cf_y()}.
}
\usage{
cos_price(K, T = 0.9, par_x, par_y, r = 0, div = 0,
          type = c("call", "put"), N = 2^10, L = 4, h = 1e-2)
}
\arguments{
  \item{K}{Numeric vector of strikes. Must be finite and strictly positive.}
  \item{T}{Time to maturity (in years).}
  \item{par_x}{List of parameters passed to \code{make_cf_x}.}
  \item{par_y}{List of parameters passed to \code{make_cf_y}.}
  \item{r}{Continuously compounded risk-free rate.}
  \item{div}{Continuously compounded dividend yield.}
  \item{type}{Option type: \code{"call"} or \code{"put"}.}
  \item{N}{Number of COS terms.}
  \item{L}{Truncation width parameter used to build the integration interval from cumulants.}
  \item{h}{Step size for finite-difference cumulant approximations from the log characteristic function.}
}
\details{
Let \eqn{X=\log S_T} with characteristic function \eqn{\varphi_X(u)}.
The function constructs
\deqn{\varphi_X(u) = \varphi_x(u)\,\varphi_y(u),}
where \eqn{\varphi_x} and \eqn{\varphi_y} are created by \code{make_cf_x(T, r, div, par_x)}
and \code{make_cf_y(T, par_y)}.

An integration interval \eqn{[a,b]} is determined from cumulants estimated by finite differences
of \eqn{\log \varphi_X(u)} using step size \code{h}. The COS expansion is then applied on \eqn{[a,b]}
to compute the option price. The function supports vector input \code{K} by evaluating each strike.
}
\value{
Numeric vector of option prices of the same length as \code{K}.
}
\references{
Fang, F. and Oosterlee, C. W. (2008).
A novel pricing method for European options based on Fourier-cosine series expansions.
\emph{SIAM Journal on Scientific Computing}, 31(2), 826--848.
}
\seealso{
\code{\link[stats]{fft}}
}
\examples{
\dontrun{
K <- c(80, 90, 100, 110, 120)
cos_price(K, T = 0.25, par_x = par_x, par_y = par_y, r = 0.01, div = 0, type = "call")
}
}
\keyword{models}
\keyword{math}
