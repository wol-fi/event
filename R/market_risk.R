cf_x <- function(u, T=1, r=0, div=0, par_x) {
  i <- 1i
  kbar <- exp(par_x$muJ + 0.5 * par_x$sigJ^2) - 1
  a <- par_x$kappa * par_x$theta
  b <- par_x$kappa - par_x$rho * par_x$sigma * i * u
  d <- sqrt(b^2 + par_x$sigma^2 * (u^2 + i * u))
  g <- (b - d) / (b + d)
  emdT <- exp(-d * T)
  C <- i * u * (r - div - par_x$lambda * kbar) * T +
    (a / par_x$sigma^2) * ((b - d) * T - 2 * log((1 - g * emdT) / (1 - g)))
  D <- ((b - d) / par_x$sigma^2) * ((1 - emdT) / (1 - g * emdT))
  jump <- exp(par_x$lambda * T * (exp(i * u * par_x$muJ - 0.5 * par_x$sigJ^2 * u^2) - 1))
  exp(i * u * log(par_x$S0x) + C + D * par_x$v0) * jump
}