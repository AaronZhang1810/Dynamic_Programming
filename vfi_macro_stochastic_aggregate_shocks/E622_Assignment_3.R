### Author: Jihuan Zhang 'Aaron' ###

# Housekeeping
rm(list = ls()); cat("\f")

# Package
library(Rtauchen) # for simulating Markov chain as approximation of AR(1) processes

# Parameters setup
alpha <- 1/3 # capital share
beta <- 143/144 # discount factor
delta <- 1/48 # depreciation
sigma <- 2 # parameter for the utility function

# Capital grid setup
kappa <- 0.4
nk <- 200 # number of capital nodes
kss <- (alpha/((1/beta)+delta-1))^(1/(1-alpha)) # steady state capital
kgrid <- seq((1-kappa)*kss, (1+kappa)*kss, length.out = nk)

# Aggregate shock grid setup
sigma_e <- 0.007 # parameter for the distribution of aggregate shocks
rho <- 0.975 # persistency of aggregate shocks
nz <- 7 # number of shock nodes
m <- 3 # max +- std. devs.
zprob <- Rtauchen(nz, sigma_e, rho, m) # simulating Markov chain that approximates the AR(1) process
zgrid <- Tgrid(nz, sigma_e, rho, m) # the transition matrix specifies transition probability from row to column
zgrid <- exp(zgrid) # since we assume log(z) ~ AR(1)

# Iteration setup
tol <- 1e-7 # tolerance
iter_max <- 1000
V0 <- matrix(0,nk,nz)
V0 <- matrix(0,nk,nz)
V1 <- matrix(0,nk,nz)
Vnew <- rep(0,nk)
C_policy <- matrix(0,nk,nz) # consumption policy function
K_policy <- matrix(0,nk,nz) # capital policy function
optimal_index <- matrix(0,nk,nz) # optimal grid position

# Value function iteration - Standard Algorithm
# this algorithm is standard and easily understandable but takes time
# for a faster algorithm, see my matlab file
for (iter in 1:iter_max) {
  for (iz in 1:nz) {
    for (ik in 1:nk) {
      for (ikp in 1:nk) {
        c <- zgrid[iz]*kgrid[ik]^alpha + (1-delta)*kgrid[ik] - kgrid[ikp]
        ifelse(c>0, Vnew[ikp] <- (c^(1-sigma)-1)/(1-sigma)+beta*sum(zprob[iz,]*V0[ikp,]), Vnew[ikp] <- -1e5)
      }
      V1[ik,iz] <- max(Vnew)
      indx <- which.max(Vnew)
      K_policy[ik,iz] = kgrid[indx]
      C_policy[ik,iz] = zgrid[iz]*kgrid[ik]^alpha+(1-delta)*kgrid[ik]-kgrid[indx]
      optimal_index[ik,iz] = ikp
    }
  }
  # we compute the vfi error by the relevant maximum distance of all
  # entries. Note that "max(abs(V1-V0))" gives a row vector, we need double max
  distance = max(max(abs((V1-V0)))/max(max(abs(V0))))
  print(iter)
  if (distance < tol) {break}
  # updating
  V0 <- V1
}

# Plotting
names <- c('z=z1','z=z2','z=z3','z=z4','z=z5','z=z6','z=z7')

matplot(as.data.frame(K_policy),type="l",xlab="k",ylab="k'(k,z)",main="Capital Policy Function")
legend("bottomright", inset=0.01, legend=names, col=c(1:7),pch=15:19,bg= ("white"), horiz=F)

matplot(as.data.frame(C_policy),type="l",xlab="k",ylab="c(k,z)",main="Consumption Policy Function")
legend("bottomright", inset=0.01, legend=names, col=c(1:7),pch=15:19,bg= ("white"), horiz=F)

matplot(as.data.frame(V1),type="l",xlab="k",ylab="V(k,z)",main="Value Function")
legend("bottomright", inset=0.01, legend=names, col=c(1:7),pch=15:19,bg= ("white"), horiz=F)