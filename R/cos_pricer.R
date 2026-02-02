cos_price <- function(K, T=0.9, r=0, div=0, par_x=list(), par_y=list(),
                      type=c("call","put"), N=1024, L=10, h=1e-4) {
  type <- match.arg(type)
  par_x <- default_x(par_x)
  par_y <- default_y(par_y)
  
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

