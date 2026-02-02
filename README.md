**Event**

## Installation

``` r
    devtools::install_github("wol-fi/event")
```

## Example
``` r
K <- 80:120
C <- cos_price(K, type="call")
plot(K, C)

Y0 <- get_Y0()
S0 <- get_S0()

```