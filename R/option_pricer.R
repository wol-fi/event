cos_price <- function(K, T=0.9, par_x, par_y, r=0, div=0,
                      type=c("call","put"), N=2^10, L=4, h=1e-2) {
  type <- match.arg(type)
  
  cfx <- make_cf_x(T, r, div, par_x)
  cfy <- make_cf_y(T, par_y)
  cf_lnS <- function(u) cfx(u) * cfy(u)
  
  cumulants <- function(cf, h=1e-4) {
    lp <- function(u) log(cf(u))
    l0 <- lp(0)
    l1p <- lp(h); l1m <- lp(-h)
    l2p <- lp(2 * h); l2m <- lp(-2 * h)
    d1 <- (l1p - l1m) / (2 * h)
    d2 <- (l1p - 2 * l0 + l1m) / (h^2)
    d4 <- (l2m - 4 * l1m + 6 * l0 - 4 * l1p + l2p) / (h^4)
    c(c1=Re(-1i * d1), c2=Re(-d2), c4=Re(d4))
  }
  
  psi_term <- function(k, c, d, a, b) {
    if (k == 0) return(d - c)
    u <- k * pi / (b - a)
    (sin(u * (d - a)) - sin(u * (c - a))) / u
  }
  
  chi_term <- function(k, c, d, a, b) {
    u <- k * pi / (b - a)
    denom <- 1 + u^2
    e1 <- cos(u * (d - a)) * exp(d) - cos(u * (c - a)) * exp(c)
    e2 <- u * (sin(u * (d - a)) * exp(d) - sin(u * (c - a)) * exp(c))
    (e1 + e2) / denom
  }
  
  vk_call <- function(k, a, b) 2/(b-a) * (chi_term(k, 0, b, a, b) - psi_term(k, 0, b, a, b))
  vk_put  <- function(k, a, b) 2/(b-a) * (psi_term(k, a, 0, a, b) - chi_term(k, a, 0, a, b))
  
  oneK <- function(Ki) {
    cs <- cumulants(cf_lnS)
    c1x <- cs["c1"] - log(Ki)
    w <- sqrt(cs["c2"] + sqrt(pmax(cs["c4"], 0)))
    a <- c1x - L * w
    b <- c1x + L * w
    if (a >= 0) a <- -1e-12
    if (b <= 0) b <-  1e-12
    
    k <- 0:(N - 1)
    u <- k * pi / (b - a)
    cfv <- cf_lnS(u)
    fk <- Re(exp(-1i * u * log(Ki)) * cfv * exp(-1i * u * a))
    fk[1] <- 0.5 * fk[1]
    vk <- if (type == "call") vapply(k, function(kk) vk_call(kk, a, b), numeric(1)) else
      vapply(k, function(kk) vk_put(kk, a, b), numeric(1))
    exp(-r * T) * Ki * sum(fk * vk)
  }
  
  sapply(K, oneK)
}

opt_price <- function(K, T=0.9, par_x, par_y, r=0, div=0,
                       type=c("call","put"), N=2^12, L=8, h=1e-4) {
  type <- match.arg(type)
  K <- as.numeric(K)
  if (any(!is.finite(K)) || any(K <= 0)) stop("K must be finite and > 0")
  if (N %% 2 != 0) stop("N must be even")
  
  cfx <- make_cf_x(T, r, div, par_x)
  cfy <- make_cf_y(T, par_y)
  cf_lnS <- function(u) cfx(u) * cfy(u)
  
  cumulants <- function(cf, h) {
    lp <- function(u) log(cf(u))
    l0 <- lp(0)
    l1p <- lp(h); l1m <- lp(-h)
    l2p <- lp(2*h); l2m <- lp(-2*h)
    d1 <- (l1p - l1m) / (2*h)
    d2 <- (l1p - 2*l0 + l1m) / (h^2)
    d4 <- (l2m - 4*l1m + 6*l0 - 4*l1p + l2p) / (h^4)
    c(c1=Re(-1i*d1), c2=Re(-d2), c4=Re(d4))
  }
  
  cs <- cumulants(cf_lnS, h=h)
  w <- sqrt(cs["c2"] + sqrt(pmax(cs["c4"], 0)))
  logK <- log(K)
  
  a0 <- cs["c1"] - L*w
  b0 <- cs["c1"] + L*w
  a <- min(a0, min(logK) - L*w)
  b <- max(b0, max(logK) + L*w)
  if (!is.finite(a) || !is.finite(b) || !(b > a)) stop("failed to build [a,b]")
  
  dx <- (b - a) / N
  du <- 2*pi / (b - a)
  
  k <- 0:(N-1)
  u <- (k - N/2) * du
  g <- cf_lnS(u) * exp(-1i * u * a)
  fhat <- fft(g)
  
  j <- 0:(N-1)
  x <- a + j * dx
  dens <- Re(fhat) * (du / (2*pi)) * ((-1)^j)
  
  dens[dens < 0] <- 0
  s0 <- sum(dens) * dx
  if (s0 > 0) dens <- dens / s0
  
  tailA <- rev(cumsum(rev(dens))) * dx
  tailB <- rev(cumsum(rev(exp(x) * dens))) * dx
  meanS <- sum(exp(x) * dens) * dx
  
  Ak <- approx(x, tailA, xout=logK, rule=2)$y
  Bk <- approx(x, tailB, xout=logK, rule=2)$y
  
  call <- exp(-r*T) * (Bk - K * Ak)
  
  if (type == "call") return(call)
  
  put <- exp(-r*T) * (K * (1 - Ak) - (meanS - Bk))
  put
}


