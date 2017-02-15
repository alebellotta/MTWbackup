set.seed(1234)

## call packages
library(stats)
library(stats4)
library(PerformanceAnalytics)
library(tseries)
## call fcts
#voaltility simulation
ornstein_uhlenbeck <- function(T.t,nstep,nsim,theta1,theta2,theta3,x0){
  dt  <- T.t/nstep
  Z <- matrix(0, ncol = nstep, nrow = nsim) #Randomnorm volatility process
  X <- matrix(0, ncol = (1+nstep), nrow = nsim) #Volatility
  X[,1] <- x0
  for (i in 2:(nstep+1)){ 
    Z[,i-1] <- rnorm(nsim, 0, sqrt(dt))
    X[,i]  <-  X[,i-1] + (theta1 - theta2*X[,i-1])*dt + theta3*sqrt(X[,i-1])*Z[,i-1]
  }
  return(OU <- list(X = X,
                    Z = Z))
}

Assetpaths <- function(S0, K, sigma, T.t, r, nsim, nstep, V, Z.V){
  dt = T.t/nstep
  S = matrix(0, nsim, (1+nstep)) 
  S[,1] = rep(S0,nsim)
  for (i in 2:(nstep+1)){ 
    eps <- rnorm(nsim)
    S[,i] <- S[,i-1]*exp((r-0.5*V[,i-1])*dt + sqrt(V[,i-1])*(rho*Z.V[,i-1] + sqrt((1-rho^2))*dt*eps))
  }
  return(S)
}

# fast algo
expBes <- function(x,nu){
  mu <- 4*nu ^2
  A1 <- 1
  A2 <- A1 * (mu - 1) / (1 * (8*x))
  A3 <- A2 * (mu - 9) / (2 * (8*x))
  A4 <- A3 * (mu - 25) / (3 * (8*x))
  A5 <- A4 * (mu - 49) / (4 * (8*x))
  A6 <- A5 * (mu - 81) / (5 * (8*x))
  A7 <- A6 * (mu -121) / (6 * (8*x))
  1/ sqrt(2*pi*x) * (A1 - A2 + A3 - A4 + A5 - A6 + A7)
}

dcCIR <- function (x, t, x0 , theta , log = FALSE ){
  c <- 2* theta[2] /((1 - exp(- theta[2] *t))* theta[3]^2)
  ncp <- 2*c*x0*exp(- theta[2] *t)     # non centrality param of chi sq
  df <- 4* theta[1] / theta[3]^2       # df of condition prob of the chi sq con. prob
  u <- c*x0* exp (- theta [2] *t)
  v <- c*x
  q <- 2* theta [1] / theta [3]^2 -1
  lik <- ( log (c) - (u+v) + q/2 * log (v/u) + log ( expBes ( 2* sqrt (u*v), q))
           + 2* sqrt (u*v))
  if(!log )
    lik <- exp(lik)
  lik
}

CIR.lik <- function(theta1 , theta2 , theta3 ) {
  n <- length(X)
  dt <- deltat(X)
  -sum(dcCIR(x=X[2: n], t=dt , x0=X[1:(n -1)] ,theta=c(theta1 , theta2 , theta3 ),
             log = TRUE ))
}

# slow algo
dcCIR2 <- function(x, t, x0 , theta , log = FALSE ){
  c <- 2* theta[2] /((1 - exp(- theta[2] *t))* theta[3]^2)
  ncp <- 2*c*x0* exp(- theta[2] *t)
  df <- 4* theta[1] / theta[3]^2
  lik <- (dchisq (2 * x * c, df = df , ncp = ncp , log = TRUE )
          + log (2*c))
  if(!log )
    lik <- exp ( lik )
  lik
}

CIR.lik2 <- function(theta1 , theta2 , theta3 ) {
  n <- length (X)
  dt <- deltat (X)
  -sum ( dcCIR2 (x=X[2: n], t=dt , x0=X[1:(n -1)] , theta =c( theta1 , theta2 , theta3 ),
                 log = TRUE ))
}



library(tseries)

AAPL = as.numeric(get.hist.quote(instrument = "KO", quote = "Close"))
ret = diff(log(AAPL))

garchpara = garch(ret, order = c(1,1))


#### stochastic variance ####
## simulate path
d <- 20                                                         # Trading days to maturity 
T.t <- d/252                                                    # Annualized time-to-maturity
nstep <- d*2                                                    # num of steps
nsim <- 10^6                                                    # num of paths
theta2 <- as.numeric(unlist(garchpara$coef[2]))                 # k
theta1 <- as.numeric(unlist(garchpara$coef[1]))*theta2          # k*long_run_mean
theta3 <- sqrt(2*theta1)-as.numeric(unlist(garchpara$coef[3]))  # volatility of the stock price volatility (2*k*theta > sigma^2)
x0 <- 0.2                                                       # sock price volatility starting point

OU <- ornstein_uhlenbeck(T.t, nstep, nsim, theta1, theta2, theta3, x0) 
V <- OU$X #sim vola matrix

# diagnostic plot
par(mfrow=c(1,2))
matplot(t(V[1:1000,]), type='l')
hist(V[,ncol(V)], breaks = 100, col='2')

## parameters calibration
X <- V[1,] #Set fit Garch(1,1)
fit <- mle(CIR.lik2, start = list(theta1 = theta1, theta2 = theta2, theta3 = theta3),
           method ="BFGS")
fit
# Coefficients:
#     theta1     theta2     theta3 
# 0.03169577 0.13281524 0.02446667  #Heston parameters

long_run_mean <- fit@coef[1]/fit@coef[2]
k <- fit@coef[2]
sigma <- fit@coef[3]
2*k*long_run_mean > sigma^2 # TRUE

OU <- ornstein_uhlenbeck(T.t, nstep, nsim, theta1 <- fit@coef[1], theta2 <- fit@coef[2], theta3 <- fit@coef[3], x0)
V <- OU$X
Z.V <- OU$Z

# diagnostic
matplot(t(V[1:10000,]), type='l')
hist(V[,ncol(V)], breaks = 100, col='2')

z <- (V[,nstep+1] - mean(V[,nstep+1]) )/sd(V[,nstep+1])
par(mfrow=c(1,2))
hist(z, breaks = 100)
hist(rnorm(nsim, 0, 1), breaks = 100)

par(mfrow=c(1,1))
hist(z, breaks = 100,  prob=T)
lines(density(rnorm(nsim, 0, 1) ), col='2')
skewness(z) ; skewness(rnorm(nsim, 0, 1))
kurtosis(z) ; kurtosis(rnorm(nsim, 0, 1))

#### log(S_t) process ####
d <- d                             # Trading days to maturity 
T.t <- T.t                         # Annualized time-to-maturity
nstep <- nstep                     # num of steps
nsim <- nsim                       # num of paths
S0 <- 25                           # S(t=0)
rho <- 0.35                        # correlation(S,V)

assetpath <- Assetpaths(S0, K, sigma, T.t, r, nsim, nstep, V, Z.V)

# diagnostic
par(mfrow=c(1,2))
matplot(t(assetpath[1:10000,]), type='l')
hist(assetpath[,ncol(assetpath)], breaks = 100, col='2')

z <- (assetpath[,nstep+1] - mean(assetpath[,nstep+1]) )/sd(assetpath[,ncol(assetpath)])
hist(z, breaks = 100)
hist(rnorm(nsim, 0, 1), breaks = 100)

par(mfrow=c(1,1))
hist(z, breaks = 100,  prob=T)
lines(density(rnorm(nsim, 0, 1) ), col='2')
skewness(z) ; skewness(rnorm(nsim, 0, 1))
kurtosis(z) ; kurtosis(rnorm(nsim, 0, 1))


#### HESTON - Version 3 ####
library(ESGtoolkit)


S0 <- tail(AAPL, n=1)                                 #Initial stock price
n <- 100000                                           #Number of simulations (feel free to reduce this)
freq <- "montly"                                     #Sampling frequency
kappa <- unlist(as.numeric(fit@coef[2]))              #Volatility mean-reversion speed
volvol <- unlist(as.numeric(fit@coef[3]))             #Volatility of volatility
rho <- -0.5                                           #Correlation between stoch. vol and spot prices
V0 <- 0.2                                            #Initial variance
theta <- unlist(as.numeric(fit@coef[1]/fit@coef[2]))  #long-term variance
r0 <- 0.0075                                           #Initial short rate

horizon <- 1# Options maturities
strikes <- 110 # Options' exercise prices

# Simulation of shocks with given correlation
set.seed(1234) # reproducibility seed
shocks <- simshocks(n =  n,
                    horizon =  horizon,
                    frequency =  freq,
                    method = "anti",
                    family = 1, par =  rho)

#  Stochastic volatility  simulation
sim.vol <- simdiff(n =  n, horizon =  horizon,
                   frequency =  freq, model = "CIR", x0 =  V0,
                   theta1 =  kappa*theta, theta2 =  kappa,
                   theta3 =  volvol, eps =  shocks[[1]])

# Stock prices simulation
sim.price <- simdiff(n = n, horizon = horizon,
                     frequency = freq, model = "GBM", x0 = S0,
                     theta1 = r0, theta2 = sqrt(sim.vol),
                     eps = shocks[[2]])


# Stock price at maturity (15 years)
S_T <- sim.price[nrow(sim.price), ]

### Monte Carlo prices
#### Estimated Monte Carlo price
discounted.payoff <- function(x)
{
  (S_T - x)*(S_T - x > 0)*exp(-r0*horizon)
}
mcprices <- sapply(strikes, function(x)
  mean(discounted.payoff(x)))

#### 95% Confidence interval around the estimation
mcprices95 <- sapply(strikes,  function(x)
  t.test(discounted.payoff(x),
         conf.level = 0.95)$conf.int)

#### 'Analytical' prices given by 'callHestoncf'
pricesAnalytic <- sapply(strikes, function(x)
  callHestoncf(S = S0, X = x, tau = horizon,
               r = r0, q = 0, v0 = V0, vT = theta,
               rho = rho, k = kappa, sigma = volvol))

results <- data.frame(cbind(strikes, mcprices,
                            t(mcprices95), pricesAnalytic))
colnames(results) <- c("strikes", "mcprices", "lower95",
                       "upper95", "pricesAnalytic")

print(results)

