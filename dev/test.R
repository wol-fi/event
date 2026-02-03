

# test: put call parity, pre-event ----------------------------------------

library(event)
rm(list=ls())

# Default: 
# - Bates model for X, see "default_x()"
# - Beta-Dist for p_T, log-normal Y_j, see "default_y()"

par_y <- default_y(list(tau=3/365, sh=0.001, sl=0.001)) 
par_x <- default_x() 

K <- seq(80, 120, 0.1)
C <- cos_price(K,T=4/365, par_x, par_y, type="call", N=2^10, L=4, h=1e-2)
P <- cos_price(K,T=4/365, par_x, par_y, type="put", N=2^10, L=4, h=1e-2)

Y0 <- get_Y0(par_y)
S0 <- get_S0(par_x, par_y)

res <- data.frame(K=K, call=C, put=P, parity_err=(C-P)-(S0-K))
plot(K, res$parity_err)

plot(K, C)
plot(K, splinefun(K, C)(K,2), type="l")


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

par_x <- default_x()
par_y <- default_y(par=list(tau=7/365))
K <- 80:120
C <- cos_price(K, T=8/365, par_x, par_y, type="call", L=6, h=1e-2)
P <- cos_price(K, T=8/365, par_x, par_y, type="put", L=6, h=1e-2)

Y0 <- get_Y0(par_y)
S0 <- get_S0(par_x, par_y)

res <- data.frame(K=K, call=C, put=P, parity_err=(C-P)-(S0-K))
plot(K, res$parity_err)




# test: pre- vs. post-event -----------------------------------------------

rm(list=ls())
library(event)

par_x <- default_x()
par_y <- default_y(par=list(tau=7/365))

Y0 <- get_Y0(par_y)
S0 <- get_S0(par_x, par_y)

K <- 80:120
C1 <- cos_price(K, T=6/365, par_x, par_y, type="call", N = 2^10, L=6, h=1e-2)
C2 <- cos_price(K, T=8/365, par_x, par_y, type="call")

plot(K, C1, type="l"); 
lines(K, C2, col=4)

plot(K, splinefun(K, C1)(K,2), type="l")
lines(K, splinefun(K, C2)(K,2), type="l", main="post-event", col=4)



# test: simulation --------------------------------------------------------

rm(list=ls())
library(event)

par_y <- default_y(par=list(tau=7/365))
sim <- wf_paths(1e3, par_y$tau, par_y$p0, par_y$sigma_p)
matplot(sim$paths[,1:10], type="l")
Yt <- get_Y_sim(sim, par_y, T=8/365)
matplot(Yt$times, Yt$dotY[,1:50], type="l")



# test: density pricer ----------------------------------------------------

rm(list=ls())
library(event)

par_x <- default_x()
par_y <- default_y(par=list(tau=7/365))

K <- 80:120
C1 <- dens_price(K, T=6/365, par_x, par_y, type="call", N = 2^10, L=6, h=1e-2)
C2 <- cos_price(K, T=8/365, par_x, par_y, type="call")

plot(K, C1, type="l"); 
lines(K, C2, col=4)

plot(K, splinefun(K, C1)(K,2), type="l")
lines(K, splinefun(K, C2)(K,2), type="l", main="post-event", col=4)

