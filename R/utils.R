

# linear payoffs ----------------------------------------------------------


get_Y0 <- function(par_y) {
  t0 <- par_y$tau
  Re(make_cf_y(t0, par_y)(-1i))
}

get_S0 <- function(par_x, par_y) {
  par_x$S0x * get_Y0(par_y)
}

# default parameter -------------------------------------------------------


default_y <- function(par=list()) {
  d <- list(
    p0=0.4,
    sigma_p=0.2,
    tau=1,
    gamma=5,
    mh=log(1.10), sh=0.02,
    ml=log(0.90), sl=0.02,
    n_quad=64
  )
  modifyList(d, par)
}

default_x <- default_x_bates <- function(par=list()) {
  d <- list(
    S0x=100,
    kappa=1.5, theta=0.04, sigma=0.4, rho=-0.7, v0=0.04,
    lambda=0.2, muJ=-0.05, sigJ=0.15
  )
  modifyList(d, par)
}

default_x_heston <- function(par=list()) {
  d <- list(
    S0x=100,
    kappa=1.5, theta=0.04, sigma=0.4, rho=-0.7, v0=0.04,
    lambda=0, muJ=0, sigJ=0
  )
  modifyList(d, par)
}

default_x_BS <- function(par=list()) {
  d <- list(
    S0x=100,
    kappa=0, theta=0, sigma=0.2, rho=-0, v0=0,
    lambda=0, muJ=0, sigJ=0
  )
  modifyList(d, par)
}


detect_model_x <- function(par_x) {
  if (length(par_x) == 0) return("bates")
  z <- function(x) isTRUE(all.equal(as.numeric(x), 0, tolerance=1e-10))
  px <- default_x_bates(par_x)
  bs <- z(px$kappa) && z(px$theta) && z(px$rho) && z(px$v0) && z(px$lambda) && z(px$muJ) && z(px$sigJ) && (abs(px$sigma) > 1e-10)
  if (bs) return("bs")
  hj <- z(px$lambda) && z(px$muJ) && z(px$sigJ)
  if (hj) return("heston")
  "bates"
}