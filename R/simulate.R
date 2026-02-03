
wf_paths <- function(n_paths, tau, p0, sigma, dt = NULL, n_steps = NULL, times = NULL, tiny = 1e-12) {
  if (!(tau > 0)) stop("tau must be > 0")
  if (!(sigma > 0)) stop("sigma must be > 0")
  if (is.null(n_steps)) {
    if (is.null(dt)) dt <- tau / 2000
    n_steps <- ceiling(tau / dt)
  }
  dt <- tau / n_steps
  if (is.null(times)) times <- seq(0, tau, by = dt)
  times <- sort(unique(pmax(0, pmin(tau, times))))
  idx <- unique(pmax(0L, pmin(n_steps, as.integer(round(times / dt)))))
  out_times <- idx * dt
  k <- length(idx)
  p <- rep(p0, length.out = n_paths)
  p <- pmax(0, pmin(1, p))
  x <- asin(2 * p - 1)
  out <- matrix(NA_real_, nrow = k, ncol = n_paths)
  pos <- 1L
  if (idx[pos] == 0L) {
    out[pos, ] <- p
    pos <- pos + 1L
  }
  sdt <- sqrt(dt)
  cap <- (pi / 2) - sqrt(tiny)
  step <- 0L
  for (j in 1:n_steps) {
    step <- step + 1L
    active <- (p > 0) & (p < 1)
    if (any(active)) {
      xa <- x[active]
      dW <- sdt * rnorm(sum(active))
      xa <- xa + 0.5 * sigma^2 * tan(xa) * dt + sigma * dW
      xa <- pmax(-cap, pmin(cap, xa))
      x[active] <- xa
      p[active] <- (1 + sin(xa)) / 2
      p[p <= tiny] <- 0
      p[p >= 1 - tiny] <- 1
      x[p == 0] <- -pi / 2
      x[p == 1] <-  pi / 2
    }
    while (pos <= k && idx[pos] == step) {
      out[pos, ] <- p
      pos <- pos + 1L
    }
  }
  list(times = out_times, paths = out)
}

get_Y_sim <- function(sim, par_y, T, tol = 1e-10, dt_post = NULL) {
  times0 <- sim$times
  p0 <- sim$paths
  if (is.null(times0) || is.null(p0)) stop("sim must have $times and $paths")
  if (!is.matrix(p0)) stop("sim$paths must be a matrix [time x path]")
  if (length(times0) != nrow(p0)) stop("length(sim$times) must equal nrow(sim$paths)")
  if (any(diff(times0) <= 0)) stop("sim$times must be strictly increasing")
  
  tau <- par_y$tau
  if (!is.finite(tau) || tau <= 0) stop("par_y$tau must be finite and > 0")
  if (!is.finite(T) || T < 0) stop("T must be finite and >= 0")
  
  tmin <- min(times0); tmax <- max(times0)
  need <- min(T, tau)
  if (need < tmin - tol || need > tmax + tol) stop("sim$times must cover min(T,tau)")
  
  gamma <- par_y$gamma
  mh <- par_y$mh; sh <- par_y$sh
  ml <- par_y$ml; sl <- par_y$sl
  
  interp_p <- function(t) {
    j_exact <- which.min(abs(times0 - t))
    if (abs(times0[j_exact] - t) <= tol) return(p0[j_exact, ])
    j <- findInterval(t, times0)
    if (j <= 0 || j >= length(times0)) stop("time not bracketed by sim$times")
    t1 <- times0[j]; t2 <- times0[j + 1]
    w <- (t - t1) / (t2 - t1)
    (1 - w) * p0[j, ] + w * p0[j + 1, ]
  }
  
  times_pre <- sort(unique(c(times0[times0 <= min(T, tau) + tol], tau)))
  if (times_pre[1] < tmin - tol) times_pre <- times_pre[times_pre >= tmin - tol]
  nP <- ncol(p0)
  p_pre <- matrix(NA_real_, nrow = length(times_pre), ncol = nP)
  for (i in seq_along(times_pre)) p_pre[i, ] <- interp_p(times_pre[i])
  
  k_h <- exp(-gamma * mh + 0.5 * gamma^2 * sh^2)
  k_l <- exp(-gamma * ml + 0.5 * gamma^2 * sl^2)
  c <- k_l / k_h
  mu_h <- exp(mh + 0.5 * (1 - 2 * gamma) * sh^2)
  mu_l <- exp(ml + 0.5 * (1 - 2 * gamma) * sl^2)
  
  qfun <- function(pp) {
    den <- pp + (1 - pp) * c
    out <- pp / den
    out[pp <= 0] <- 0
    out[pp >= 1] <- 1
    out
  }
  
  ybar <- function(pp) {
    qq <- qfun(pp)
    mu_l + (mu_h - mu_l) * qq
  }
  
  dotY_pre <- ybar(p_pre)
  
  if (T <= tau + tol) {
    keep <- which(times_pre <= T + tol)
    return(list(times = times_pre[keep], dotY = dotY_pre[keep, , drop = FALSE]))
  }
  
  itau <- which.min(abs(times_pre - tau))
  p_tau <- p_pre[itau, ]
  is_h <- runif(nP) < p_tau
  lnY <- numeric(nP)
  if (any(is_h))  lnY[is_h]  <- rnorm(sum(is_h),  mean = mh, sd = sh)
  if (any(!is_h)) lnY[!is_h] <- rnorm(sum(!is_h), mean = ml, sd = sl)
  Y_tau <- exp(lnY)
  
  if (is.null(dt_post)) {
    dts <- diff(times0)
    dt_post <- stats::median(dts[dts > 0])
  }
  times_post <- seq(from = tau, to = T, by = dt_post)
  if (abs(tail(times_post, 1) - T) > tol) times_post <- c(times_post, T)
  times_out <- sort(unique(c(times_pre[times_pre < tau - tol], tau, times_post)))
  
  dotY <- matrix(NA_real_, nrow = length(times_out), ncol = nP)
  pre_idx <- which(times_out < tau - tol)
  if (length(pre_idx) > 0) {
    map <- match(times_out[pre_idx], times_pre)
    dotY[pre_idx, ] <- dotY_pre[map, , drop = FALSE]
  }
  post_idx <- which(times_out >= tau - tol)
  dotY[post_idx, ] <- matrix(rep(Y_tau, each = length(post_idx)), nrow = length(post_idx))
  
  list(times = times_out, dotY = dotY, Y_tau = Y_tau, is_h = is_h)
}
