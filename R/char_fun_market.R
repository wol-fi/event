cf_x_bates <- function(u, par) {
  i <- 1i
  kbar <- exp(par$muJ + 0.5 * par$sigJ^2) - 1
  a <- par$kappa * par$theta
  b <- par$kappa - par$rho * par$sigma * i * u
  d <- sqrt(b^2 + par$sigma^2 * (u^2 + i * u))
  g <- (b - d) / (b + d)
  emdT <- exp(-d * par$T)
  C <- i * u * (par$r - par$div - par$lambda * kbar) * par$T +
    (a / par$sigma^2) * ((b - d) * par$T - 2 * log((1 - g * emdT) / (1 - g)))
  D <- ((b - d) / par$sigma^2) * ((1 - emdT) / (1 - g * emdT))
  jump <- exp(par$lambda * par$T * (exp(i * u * par$muJ - 0.5 * par$sigJ^2 * u^2) - 1))
  exp(i * u * log(par$S0x) + C + D * par$v0) * jump
}