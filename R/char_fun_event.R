cf_y <- function(u, par) {
  p0 <- par$p0
  vp <- par$p_var
  if (p0 <= 0 || p0 >= 1) stop("p0 must be in (0,1)")
  vmax <- p0 * (1 - p0)
  if (vp <= 0 || vp >= vmax) stop("p_var must satisfy 0 < p_var < p0*(1-p0)")
  
  gamma <- par$gamma
  mh <- par$mh; ml <- par$ml
  sh <- par$sh; sl <- par$sl
  
  s <- p0 * (1 - p0) / vp - 1
  a <- p0 * s
  b <- (1 - p0) * s
  
  alpha <- b - 1
  beta <- a - 1
  n <- par$n_quad
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
  eg <- eigen(J, symmetric = TRUE)
  pk <- (eg$values + 1) / 2
  wk <- eg$vectors[1, ]^2
  wk <- wk / sum(wk)
  
  k_h <- exp(-gamma * mh + 0.5 * (gamma^2) * sh^2)
  k_l <- exp(-gamma * ml + 0.5 * (gamma^2) * sl^2)
  mu_h <- exp(mh + 0.5 * (1 - 2 * gamma) * sh^2)
  mu_l <- exp(ml + 0.5 * (1 - 2 * gamma) * sl^2)
  
  c <- k_l / k_h
  
  wk2 <- wk * (k_l + (k_h - k_l) * pk)
  denom <- sum(wk2)
  
  qk <- pk / (pk + (1 - pk) * c)
  xk <- mu_l + (mu_h - mu_l) * qk
  
  if (length(u) == 1) sum(wk2 * exp(1i * u * log(xk))) / denom else
    sapply(u, function(uu) sum(wk2 * exp(1i * uu * log(xk))) / denom)
}
