\name{wf_paths}
\alias{wf_paths}
\title{Simulate Wright--Fisher probability paths on [0,1]}
\description{
Simulates sample paths of a Wright--Fisher-type diffusion for an outcome probability process
\eqn{p_t \in [0,1]} on \eqn{[0,\tau]}. The implementation evolves a transformed state
\eqn{x_t = \arcsin(2p_t-1)} with an Euler--Maruyama step and clamps near the boundaries to
avoid numerical blow-ups.
}
\usage{
wf_paths(n_paths, tau, p0, sigma, dt = NULL, n_steps = NULL, times = NULL, tiny = 1e-12)
}
\arguments{
  \item{n_paths}{Integer number of Monte Carlo paths.}
  \item{tau}{Positive scalar horizon \eqn{\tau}.}
  \item{p0}{Initial probability \eqn{p_0}. Values are clamped to \eqn{[0,1]}.}
  \item{sigma}{Positive scalar diffusion scale parameter.}
  \item{dt}{Optional time step. If \code{n_steps} is \code{NULL}, defaults to \code{tau/2000}.}
  \item{n_steps}{Optional integer number of steps on \eqn{[0,\tau]}. If provided, \code{dt} is set to \code{tau/n_steps}.}
  \item{times}{Optional numeric vector of requested output times in \eqn{[0,\tau]}. If \code{NULL}, outputs the full simulation grid.}
  \item{tiny}{Small positive value used to snap probabilities close to 0 or 1 exactly to 0 or 1, and to define boundary caps in the transform.}
}
\details{
The function simulates on a uniform grid of \code{n_steps} with \code{dt = tau/n_steps}.
Requested output times are projected onto this grid and returned as \code{$times}.

Internally it uses the transform \eqn{x=\arcsin(2p-1)} and updates, for active paths with
\eqn{p\in(0,1)}, via an Euler step of the form
\deqn{x_{t+\Delta t} \approx x_t + \tfrac12\sigma^2 \tan(x_t)\Delta t + \sigma \Delta W_t,}
then maps back by \eqn{p=(1+\sin x)/2}. Paths are clamped to a cap strictly inside
\eqn{(-\pi/2,\pi/2)} determined by \code{tiny}, and probabilities within \code{tiny} of 0 or 1 are snapped
to the boundary.
}
\value{
A list with components:
\itemize{
  \item \code{times}: numeric vector of output times.
  \item \code{paths}: numeric matrix of dimension \code{length(times) x n_paths} containing simulated \eqn{p_t} values.
}
}
\note{
This is a discretization-based simulation. For large \code{sigma} or coarse \code{dt}, boundary behavior and moments can be sensitive to discretization and clamping choices.
}
\seealso{
\code{\link{get_Y_sim}}
}
\examples{
\dontrun{
sim <- wf_paths(n_paths = 1000, tau = 1/12, p0 = 0.4, sigma = 1.0, n_steps = 2000)
str(sim)
}
}
\keyword{models}
\keyword{math}