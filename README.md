**Event**

## Installation

``` r
    devtools::install_github("wol-fi/event")
```

## Examples

### A) Pricing European Call: pre- and post-event expiry
- Market risk `X`: Bates model (stochastic vol + Merton jumps)
- Event risk `Y`: stochastic p_t (Beta dist.) + mix of log-normals
- check parameters (`par_x`, `par_y`) for details

<img width="800" height="400" alt="call" src="https://github.com/user-attachments/assets/ab2f8dbf-44a0-4d04-b214-5830bf7dbf10" />
<img width="800" height="400" alt="rnd" src="https://github.com/user-attachments/assets/7281b14b-967e-4800-bfe3-eb10ced42000" />
<img width="800" height="400" alt="iv" src="https://github.com/user-attachments/assets/c65797e4-d447-42d4-9fa3-f933a54c4f00" />

``` r
library(event)

par_x <- default_x()
par_y <- default_y(par=list(tau=7/365))

Y0 <- get_Y0(par_y)
S0 <- get_S0(par_x, par_y)

K <- 80:120
C1 <- opt_price(K, T=6/365, par_x, par_y, type="call") # option price via RND
C2 <- cos_price(K, T=8/365, par_x, par_y, type="call") # option price via COS-method

plot(K, C1, type="l", main="Call price") 
lines(K, C2, col=4)
legend("topright", legend=c("pre-event", "post-event"), lty=1, col=c(1,4), bty="n")

plot(K, splinefun(K, C1)(K,2), type="l", ylab="RND", main="RND")
lines(K, splinefun(K, C2)(K,2), type="l", col=4)
legend("topright", legend=c("pre-event", "post-event"), lty=1, col=c(1,4), bty="n")

iv1 <- ivol(K, C1, S0, 6/365, type="call")
iv2 <- ivol(K, C2, S0, 8/365, type="call")
plot(log(K/S0), iv1, type="l", ylim=range(c(iv1, iv2)), ylab="iv", main="implied vola")
lines(log(K/S0), iv2, col=4)
legend("topright", legend=c("pre-event", "post-event"), lty=1, col=c(1,4), bty="n")
```

### B) Simulation of Event-Risk Multiplier

<img width="800" height="400" alt="Yt" src="https://github.com/user-attachments/assets/3f55961d-8048-4847-8758-a0889e8b3c26" />

``` r
library(event)

par_y <- default_y(par=list(tau=7/365, sigma_p=2))
sim <- wf_paths(100, par_y$tau, par_y$p0, par_y$sigma_p)
matplot(sim$paths, type="l")
Yt <- get_Y_sim(sim, par_y, T=8/365)
matplot(Yt$times*365, Yt$dotY, type="l", main="Event-Risk Multiplier", xlab="days", ylab="Y_t")

``` 

### C) Impact of `sigma_p`

<img width="800" height="400" alt="sigma_p" src="https://github.com/user-attachments/assets/fc448b45-4a8a-4e73-890b-6930daeb1ef9" />

``` r
library(event)

par_x <- default_x()
par_y1 <- par_y2 <- default_y(list(tau=1, sh=0.1, sl=0.1, p0=0.5, n_quad=2^8))
par_y1$sigma_p <- 0.6
par_y2$sigma_p <- 0.001

K <- 80:140
C1 <- opt_price(K, T=0.95, par_x, par_y1, type="call")
C2 <- opt_price(K, T=0.95, par_x, par_y2, type="call")
S0 <- get_S0(par_x, par_y1)

iv1 <- ivol(K, C1, S0, 0.9)
iv2 <- ivol(K, C2, S0, 0.9)
plot(log(K/S0), iv1, type="l", ylim=range(c(iv1, iv2)), main="implied vola", ylab="iv"); grid()
lines(log(K/S0), iv2, col=4)
legend("topright", legend=c(paste0("sigma_p=",par_y1$sigma_p), paste0("sigma_p=",par_y2$sigma_p )), lty=1, col=c(1,4), bty="n")

```
