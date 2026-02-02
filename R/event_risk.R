make_cf_y <- function(T, par_y=list()) {
  par_y <- default_y(par_y)
  
  p0 <- par_y$p0
  sigma_p <- par_y$sigma_p
  tau <- par_y$tau
  gamma <- par_y$gamma
  mh <- par_y$mh; sh <- par_y$sh
  ml <- par_y$ml; sl <- par_y$sl
  n <- par_y$n_quad
  
  if (!(p0 > 0 && p0 < 1)) stop("par_y$p0 must be in (0,1)")
  if (!(T >= 0)) stop("T must be >= 0")
  if (!(tau > 0)) stop("par_y$tau must be > 0")
  
  k_h <- exp(-gamma * mh + 0.5 * gamma^2 * sh^2)
  k_l <- exp(-gamma * ml + 0.5 * gamma^2 * sl^2)
  c <- k_l / k_h
  
  Eh <- function(u) exp((1i * u - gamma) * mh + 0.5 * (1i * u - gamma)^2 * sh^2)
  El <- function(u) exp((1i * u - gamma) * ml + 0.5 * (1i * u - gamma)^2 * sl^2)
  
  if (T >= tau) {
    q0 <- p0 / (p0 + (1 - p0) * c)
    f <- function(u) q0 * Eh(u) / k_h + (1 - q0) * El(u) / k_l
    return(function(u) if (length(u) == 1) f(u) else sapply(u, f))
  }
  
  if (!(sigma_p > 0 && T > 0)) {
    pk <- p0
    wk <- 1
  } else {
    den <- exp(sigma_p^2 * T) - 1
    a_par <- p0 / den
    b_par <- (1 - p0) / den
    alpha <- b_par - 1
    beta <- a_par - 1
    
    k <- 1:n
    A <- (beta^2 - alpha^2) / ((2 * k + alpha + beta) * (2 * k + alpha + beta + 2))
    B <- sqrt(4 * k * (k + alpha) * (k + beta) * (k + alpha + beta) /
                ((2 * k + alpha + beta)^2 * (2 * k + alpha + beta + 1) * (2 * k + alpha + beta - 1)))
    J <- matrix(0, n, n)
    diag(J) <- A
    if (n > 1) {
      diag(J[-1, -n]) <- B[-n]
      diag(J[-n, -1]) <- B[-n]
    }
    eg <- eigen(J, symmetric=TRUE)
    pk <- (eg$values + 1) / 2
    wk <- eg$vectors[1, ]^2
    wk <- wk / sum(wk)
  }
  
  wk2 <- wk * (k_l + (k_h - k_l) * pk)
  denom <- sum(wk2)
  qk <- pk / (pk + (1 - pk) * c)
  
  mu_h <- exp(mh + 0.5 * (1 - 2 * gamma) * sh^2)
  mu_l <- exp(ml + 0.5 * (1 - 2 * gamma) * sl^2)
  xk <- mu_l + (mu_h - mu_l) * qk
  lxk <- log(xk)
  
  function(u) {
    if (length(u) == 1) sum(wk2 * exp(1i * u * lxk)) / denom else
      sapply(u, function(uu) sum(wk2 * exp(1i * uu * lxk)) / denom)
  }
}

cf_y <- function(u, T, par_y=list()) make_cf_y(T, par_y)(u)

# -------------------------------------------------------------------------

get_Y0 <- function(par_y=list()) {
  par_y <- default_Y(par_y)
  t0 <- par_y$tau
  Re(make_cf_y(t0, par_y)(-1i))
}

get_S0 <- function(T, par_y=list(), par_x=list()) {
  par_x$S0x * get_Y0(T, par_y)
}