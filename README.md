**Event**

## Installation

``` r
    devtools::install_github("wol-fi/event")
```

## Example
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

plot(K, splinefun(K, C1)(K,2), type="l", ylab="RND", main="RND")
lines(K, splinefun(K, C2)(K,2), type="l", col=4)

```