

# test: put call parity ----------------------------------------------------

library(event)
rm(list=ls())

# Default: 
# - Bates model for X, see "default_x()"
# - Beta-Dist for p_T, log-normal Y_j, see "default_y()"

K <- 80:120
C <- cos_price(K, type="call")
P <- cos_price(K, type="put")

Y0 <- get_Y0()
S0 <- get_S0()

res <- data.frame(K=K, call=C, put=P, parity_err=(C-P)-(S0-K))
plot(K, res$parity_err)



# test: different p_var ---------------------------------------------------

library(event)
rm(list=ls())

par_y1 <- par_y2 <- default_y()
par_y2$p_var <- 0.001

K <- 80:120
C1 <- cos_price(K, par_y=par_y1, type="call")
C2 <- cos_price(K, par_y=par_y2, type="call")

plot(K, C1, type="l"); grid()
lines(K, C2, col=4)
