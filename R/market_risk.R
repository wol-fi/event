make_cf_x <- function(T=1, r=0, div=0, par_x) {
  model <- detect_model_x(par_x)
  
  S0x <- par_x$S0x
  kappa <- par_x$kappa
  theta <- par_x$theta
  sigma <- par_x$sigma
  rho <- par_x$rho
  v0 <- par_x$v0
  lambda <- par_x$lambda
  muJ <- par_x$muJ
  sigJ <- par_x$sigJ
  
  ls0 <- log(S0x)
  
  if (model == "bs") {
    vol <- sigma
    function(u) {
      i <- 1i
      exp(i*u*(ls0 + (r - div - 0.5*vol^2)*T) - 0.5*vol^2*u^2*T)
    }
  } else {
    kbar <- exp(muJ + 0.5*sigJ^2) - 1
    a <- kappa * theta
    function(u) {
      i <- 1i
      b <- kappa - rho * sigma * i * u
      d <- sqrt(b^2 + sigma^2 * (u^2 + i * u))
      g <- (b - d) / (b + d)
      emdT <- exp(-d * T)
      C <- i * u * (r - div - lambda * kbar) * T +
        (a / sigma^2) * ((b - d) * T - 2 * log((1 - g * emdT) / (1 - g)))
      D <- ((b - d) / sigma^2) * ((1 - emdT) / (1 - g * emdT))
      if (model == "heston") {
        exp(i * u * ls0 + C + D * v0)
      } else {
        jump <- exp(lambda * T * (exp(i * u * muJ - 0.5 * sigJ^2 * u^2) - 1))
        exp(i * u * ls0 + C + D * v0) * jump
      }
    }
  }
}

cf_x <- function(u, T=1, r=0, div=0, par_x) {
  make_cf_x(T, r, div, par_x)(u)
}

