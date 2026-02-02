


# test: put call parity, pre-event ----------------------------------------

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

plot(K, C)

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


# test: put call parity, post-event ---------------------------------------


library(event)
rm(list=ls())

K <- 60:120
C <- cos_price(K, T=.9/365, tau=1/365, type="call")
P <- cos_price(K, T=.9/365, tau=1/365, type="put")

Y0 <- get_Y0()
S0 <- get_S0()

res <- data.frame(K=K, call=C, put=P, parity_err=(C-P)-(S0-K))
plot(K, res$parity_err)




# test: pre- vs. post-event -----------------------------------------------

library(event)
rm(list=ls())

par_y <- default_y()
par_x <- default_x()
# par_y$sl <- 0.0001
# par_y$sh <- 0.0001
par_y$tau <- 7/365
K <- 80:120
C1 <- cos_price(K, T=6/365, par_y=par_y, type="call", N = 2^10, L=6, h=1e-2)
C2 <- cos_price(K, T=8/365, par_y=par_y, type="call")

plot(K, C1, type="l"); 
lines(K, C2, col=4)

plot(K, splinefun(K, C1)(K,2), type="l")
lines(K, splinefun(K, C2)(K,2), type="l", main="post-event", col=4)

