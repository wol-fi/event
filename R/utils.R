default_y <- function(par=list()) {
  d <- list(
    p0=0.4,
    p_var=NA_real_,
    gamma=5,
    mh=log(1.10), sh=0.10,
    ml=log(0.90), sl=0.10,
    n_quad=64
  )
  par <- modifyList(d, par)
  if (!is.finite(par$p_var)) par$p_var <- 0.6 * par$p0 * (1 - par$p0)
  par
}

default_x <- function(par=list()) {
  d <- list(
    S0x=100,
    kappa=1.5, theta=0.04, sigma=0.4, rho=-0.7, v0=0.04,
    lambda=0.2, muJ=-0.05, sigJ=0.15
  )
  modifyList(d, par)
}