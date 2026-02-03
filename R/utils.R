

# linear payoffs ----------------------------------------------------------


get_Y0 <- function(par_y) {
  t0 <- par_y$tau
  Re(make_cf_y(t0, par_y)(-1i))
}

get_S0 <- function(par_x, par_y) {
  par_x$S0x * get_Y0(par_y)
}


# ivol --------------------------------------------------------------------


ivol <- function(K, cp_price, S0, T, r=0, div=0, type=c("call","put"),
                 tol=1e-10, max_iter=200, v_max=10, v0=0.2) {
  type <- match.arg(type)
  K <- as.numeric(K)
  p <- as.numeric(cp_price)
  if (length(p) == 1L) p <- rep(p, length(K))
  if (length(K) != length(p)) stop("length(K) must equal length(cp_price) (or cp_price scalar)")
  if (!is.finite(S0) || S0 <= 0) stop("S0 must be finite and > 0")
  if (!is.finite(T) || T <= 0) stop("T must be finite and > 0")
  if (!is.finite(r) || !is.finite(div)) stop("r and div must be finite")
  if (any(!is.finite(K)) || any(K <= 0)) stop("K must be finite and > 0")
  
  df_r <- exp(-r*T)
  F0 <- S0 * exp((r - div)*T)
  
  bs_call_fwd <- function(x, v) {
    if (!is.finite(v) || v <= 0) return(pmax(1 - exp(x), 0))
    pnorm(-x/v + v/2) - exp(x) * pnorm(-x/v - v/2)
  }
  
  call_from_put <- function(P, K) {
    C <- P + df_r * (F0 - K)
    C
  }
  
  out <- rep(NA_real_, length(K))
  
  for (i in seq_along(K)) {
    Ki <- K[i]
    if (!is.finite(p[i])) next
    
    x <- log(Ki / F0)
    
    C_disc <- if (type == "call") p[i] else call_from_put(p[i], Ki)
    c_fwd <- C_disc / (df_r * F0)
    
    c0 <- pmax(1 - exp(x), 0)
    c1 <- 1
    if (!(is.finite(c_fwd) && c_fwd >= c0 - 10*tol && c_fwd <= c1 + 10*tol)) next
    
    if (c_fwd <= c0 + tol) {
      out[i] <- 0
      next
    }
    
    lo <- 0
    hi <- v0
    phi <- bs_call_fwd(x, hi)
    while (is.finite(phi) && phi < c_fwd && hi < v_max) {
      hi <- 2 * hi
      phi <- bs_call_fwd(x, hi)
    }
    if (!is.finite(phi) || phi < c_fwd) next
    
    for (iter in 1:max_iter) {
      mid <- 0.5 * (lo + hi)
      pmid <- bs_call_fwd(x, mid)
      if (!is.finite(pmid)) { lo <- NA_real_; break }
      if (abs(pmid - c_fwd) <= tol) { lo <- mid; hi <- mid; break }
      if (pmid < c_fwd) lo <- mid else hi <- mid
      if (hi - lo <= tol * pmax(1, mid)) break
    }
    
    v_imp <- 0.5 * (lo + hi)
    out[i] <- v_imp / sqrt(T)
  }
  
  out
}



# default parameter -------------------------------------------------------


default_y <- function(par=list()) {
  d <- list(
    p0=0.4,
    sigma_p=0.2,
    tau=7/365,
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