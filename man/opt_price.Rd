\name{opt_price}
\alias{opt_price}
\title{Density-based European option pricing from a characteristic function}
\description{
Prices European call or put options for a vector of strikes by first inverting the
characteristic function of the log-price to a log-density on a uniform grid (FFT-based),
then using precomputed tail integrals to evaluate option prices efficiently across strikes.
}
\usage{
opt_price(K, T = 0.9, par_x, par_y, r = 0, div = 0,
           type = c("call", "put"), N = 2^12, L = 8, h = 1e-4)
}
\arguments{
  \item{K}{Numeric vector of strikes. Must be finite and strictly positive.}
  \item{T}{Time to maturity (in years).}
  \item{par_x}{List of parameters passed to \code{make_cf_x}.}
  \item{par_y}{List of parameters passed to \code{make_cf_y}.}
  \item{r}{Continuously compounded risk-free rate.}
  \item{div}{Continuously compounded dividend yield.}
  \item{type}{Option type: \code{"call"} or \code{"put"}.}
  \item{N}{Number of grid points for the log-density inversion. Must be even.}
  \item{L}{Truncation width parameter used to form the log-price interval \eqn{[a,b]} from cumulants.}
  \item{h}{Step size for finite-difference cumulant approximations from the log characteristic function.}
}
\details{
Let \eqn{X=\log S_T} with characteristic function \eqn{\varphi_X(u)}. The function builds an interval
\eqn{[a,b]} from cumulants estimated via finite differences of \eqn{\log \varphi_X(u)}, then inverts
\eqn{\varphi_X} to an approximation of the density \eqn{f_X} on a uniform grid by an FFT-consistent
Fourier inversion.

Given \eqn{k=\log K}, call prices are computed as
\deqn{C(K) = e^{-rT}\left(\int_k^\infty e^x f_X(x)\,dx - K \int_k^\infty f_X(x)\,dx\right).}
The two tail integrals are precomputed on the grid via reverse cumulative sums and interpolated to the
requested log-strikes. Put prices are obtained from the same tail integrals using
\deqn{P(K)=e^{-rT}\left(K(1-A(k))-(\mathbb{E}[e^X]-B(k))\right),}
where \eqn{A(k)=\int_k^\infty f_X(x)\,dx} and \eqn{B(k)=\int_k^\infty e^x f_X(x)\,dx}.
}
\value{
Numeric vector of option prices of the same length as \code{K}.
}
\references{
Carr, P. and Madan, D. B. (1999).
Option valuation using the fast Fourier transform.
\emph{Journal of Computational Finance}, 2(4), 61--73.

Fang, F. and Oosterlee, C. W. (2008).
A novel pricing method for European options based on Fourier-cosine series expansions.
\emph{SIAM Journal on Scientific Computing}, 31(2), 826--848.
}
\seealso{
\code{\link{cos_price}} (if provided in the same package),
\code{\link[stats]{fft}},
\code{\link[stats]{approx}}.
}
\examples{
\dontrun{
# Requires make_cf_x() and make_cf_y() to be defined.
# K <- c(80, 90, 100, 110, 120)
# opt_price(K, T = 0.25, par_x = par_x, par_y = par_y, r = 0.01, div = 0, type = "call")
}}
\keyword{models}
\keyword{math}
