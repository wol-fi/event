**Event**

## Installation

``` r
    devtools::install_github("wol-fi/event")
```

## Example
- Pricing European call options with expiry before- and after event
- Market risk (X): Bates model (stochastic vol + Merton jumps)
- Event risk (Y): stochastic p_t (Beta dist.) + mix of log-normals
- check parameters (par_x, par_y) for details

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
```
<img width="800" height="400" alt="call" src="https://github.com/user-attachments/assets/ab2f8dbf-44a0-4d04-b214-5830bf7dbf10" />

``` r
plot(K, splinefun(K, C1)(K,2), type="l", ylab="RND", main="RND")
lines(K, splinefun(K, C2)(K,2), type="l", col=4)
legend("topright", legend=c("pre-event", "post-event"), lty=1, col=c(1,4), bty="n")
```
<img width="800" height="400" alt="rnd" src="https://github.com/user-attachments/assets/7281b14b-967e-4800-bfe3-eb10ced42000" />

``` r
iv1 <- ivol(K, C1, S0, 6/365, type="call")
iv2 <- ivol(K, C2, S0, 8/365, type="call")
plot(log(K/S0), iv1, type="l", ylim=range(c(iv1, iv2)), ylab="iv", main="implied vola")
lines(log(K/S0), iv2, col=4)
legend("topright", legend=c("pre-event", "post-event"), lty=1, col=c(1,4), bty="n")
```
<img width="800" height="400" alt="iv" src="https://github.com/user-attachments/assets/c65797e4-d447-42d4-9fa3-f933a54c4f00" />


