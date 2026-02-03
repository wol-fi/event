\name{get_Y_sim}
\alias{get_Y_sim}
\title{Compute simulated event-risk multipliers from probability paths}
\description{
Given simulated probability paths \eqn{p_t} (typically from \code{\link{wf_paths}}), computes the
simulated event-risk multiplier process \eqn{\dot Y_t} up to maturity \eqn{T}, with a single event time
\eqn{\tau} supplied in \code{par_y}. For \eqn{t<\tau}, \eqn{\dot Y_t} is the conditional expectation of the
event multiplier under the tilted measure induced by \code{par_y}. For \eqn{t\ge \tau} the multiplier is
constant along each path and equal to the realized event draw at \eqn{\tau}.
}
\usage{
get_Y_sim(sim, par_y, T, tol = 1e-10, dt_post = NULL)
}
\arguments{
  \item{sim}{List with \code{$times} (strictly increasing numeric vector) and \code{$paths} (matrix \code{[time x path]}) containing simulated \eqn{p_t}.}
  \item{par_y}{List of event parameters. Must contain \code{tau}, \code{gamma}, \code{mh}, \code{sh}, \code{ml}, \code{sl}.}
  \item{T}{Maturity time \eqn{T\ge 0}. The function requires that \code{sim$times} covers \eqn{\min(T,\tau)}.}
  \item{tol}{Numerical tolerance used for matching/interpolation and time comparisons.}
  \item{dt_post}{Optional post-event output step size for \eqn{t\in[\tau,T]}. If \code{NULL}, uses the median positive increment of \code{sim$times}.}
}
\details{
Let \eqn{p_t} denote the simulated probability (under the simulation measure used for \code{sim}).
Given \code{par_y}, define
\deqn{k_h = \exp(-\gamma m_h + \tfrac12\gamma^2 s_h^2),\quad k_l = \exp(-\gamma m_l + \tfrac12\gamma^2 s_l^2),\quad c = k_l/k_h.}
The function forms a tilted probability
\deqn{q(p) = \frac{p}{p + (1-p)c},}
clipped to \eqn{[0,1]} at the boundaries. It then computes the pre-event conditional mean multiplier
\deqn{\dot Y_t = \bar y(p_t) = \mu_l + (\mu_h-\mu_l)\,q(p_t),}
with
\deqn{\mu_h = \exp(m_h + \tfrac12(1-2\gamma)s_h^2),\qquad \mu_l = \exp(m_l + \tfrac12(1-2\gamma)s_l^2).}

If \eqn{T\le \tau}, only the pre-event values up to \eqn{T} are returned.

If \eqn{T>\tau}, the function draws a Bernoulli event indicator at \eqn{\tau} with success probability \eqn{p_\tau} (pathwise),
samples a lognormal shock conditional on the outcome using \code{mh,sh} or \code{ml,sl},
sets \eqn{Y_\tau=\exp(\text{shock})}, and outputs \eqn{\dot Y_t=Y_\tau} for all \eqn{t\in[\tau,T]} along each path.
}
\value{
A list with components:
\itemize{
  \item \code{times}: numeric vector of output times (includes \eqn{\tau} and, if \eqn{T>\tau}, a post-event grid up to \eqn{T}).
  \item \code{dotY}: numeric matrix \code{[length(times) x n_paths]} of simulated \eqn{\dot Y_t}.
  \item \code{Y_tau}: (only if \eqn{T>\tau}) numeric vector of realized multipliers \eqn{Y_\tau} per path.
  \item \code{is_h}: (only if \eqn{T>\tau}) logical vector indicating high-state draws at \eqn{\tau}.
}
}
\seealso{
\code{\link{wf_paths}}
}
\examples{
\dontrun{
sim <- wf_paths(n_paths = 1000, tau = 1/12, p0 = 0.4, sigma = 1.0, n_steps = 2000)
par_y <- list(tau = 1/12, gamma = 5, mh = log(1.1), sh = 0.1, ml = log(0.9), sl = 0.1)
ys <- get_Y_sim(sim, par_y, T = 0.2)
str(ys)
}
}
\keyword{models}
\keyword{math}