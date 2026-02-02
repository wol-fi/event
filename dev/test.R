

# test put call parity ----------------------------------------------------

library(event)

p0 <- 0.4
p_var <- 0.6 * p0*(1-p0) 
par1 <- list(T=1,S0x=100,r=0,div=0,kappa=1.5,theta=0.04,sigma=0.4,rho=-0.7,v0=0.04,lambda=0.2,muJ=-0.05,sigJ=0.15,
             p0=p0, p_var=p_var,gamma=5,mh=log(1.10),sh=0.10,ml=log(0.90),sl=0.10,N=1024,L=10,n_quad=64)

K <- 80:120
C <- cos_price(K, par1, "call")
