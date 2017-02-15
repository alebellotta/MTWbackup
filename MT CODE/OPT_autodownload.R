####  DATA  ####
#get risk free rate#
getSymbols.FRED("TB1YR", env = globalenv())
TB1YR = as.data.frame(TB1YR)
TB1YR = TB1YR$TB1YR
r = tail(TB1YR, n=1)

#get data#
ticker = "AAPL"
stockprice = as.numeric(get.hist.quote(instrument = ticker, quote = "Close"))


start.time <- Sys.time()
ret = diff(log(stockprice))
S0 = tail(stockprice, n=1) 


K = S0*1.2
T = 2
tau = T

sigma = sd(ret[length(ret)-250:length(ret)])
n = 1000
eps = rnorm(n, mean = 0, sd = sigma)
Nsim = 1000
Nsteps = 1000
Npilot = 1000
Sb = K

Call = GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, b = r, sigma = sigma)

############################
####  Closed Formulas   ####
############################

BScall=function(S0,K,sigma,T,r){
  d1=(log(S0/K)+(r+0.5*sigma^2)*T)/(sigma*sqrt(T))
  d2=d1-sigma*sqrt(T);
  Price=S0*pnorm(d1)-K*exp(-r*T)*pnorm(d2)
  return(Price)
}

######################
####   Crude MC   ####
######################

BScallMC=function(S0,K,sigma,T,r,n){
  S=S0*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*rnorm(n))
  disc_payoff=exp(-r*T)*(ifelse(S>K,S-K,0))
  return(list(mean(disc_payoff),sd(disc_payoff)))
}

################################
####  Antithetic  variates  ####
################################

BScallMCav=function(S0,K,sigma,T,r,n){
  eps=rnorm(n)
  S1=S0*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*eps)
  S2=S0*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*(-eps))
  payoff1=ifelse(S1>K,S1-K,0)
  payoff2=ifelse(S2>K,S2-K,0)
  disc_payoff=exp(-r*T)*0.5*(payoff1+payoff2)
  return(list(mean(disc_payoff),sd(disc_payoff)))
}

############################
####  Control variates  ####
############################

BScallMCcv=function(S0,K,sigma,T,r,n,npilot){
  Spilot=S0*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*rnorm(Npilot))
  Cpilot=exp(-r*T)*ifelse(Spilot>K,Spilot-K,0)
  Cov=cov(Spilot,Cpilot)
  VarS=S0^2*exp(2*r*T)*(exp(T*sigma^2)-1)
  c=-Cov/VarS
  ExpS=S0*exp(r*T)
  S=S0*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*rnorm(n))
  C=exp(-r*T)*ifelse(S>K,S-K,0)
  controlvar=C+c*(S-ExpS)
  return(list(mean(controlvar),sd(controlvar)))
}

##################
####  Results ####
##################
Cbs=BScall(S0,K,sigma,T,r)
Cbsmc=BScallMC(S0,K,sigma,T,r,n)
Cbsmcav=BScallMCav(S0,K,sigma,T,r,n)
Cbsmccv=BScallMCcv(S0,K,sigma,T,r,n,npilot)


####  HESTON  ####
set.seed(1234)

## call fcts
#volatility simulation
ornstein_uhlenbeck = function(T.t,nstep,nsim,theta1,theta2,theta3,x0){
  dt = T.t/nstep
  Z = matrix(0, ncol = nstep, nrow = nsim) #Randomnorm volatility process
  X = matrix(0, ncol = (1+nstep), nrow = nsim) #Volatility
  X[,1] = x0
  for (i in 2:(nstep+1)){ 
    Z[,i-1] = rnorm(nsim, 0, sqrt(dt))
    X[,i]  =  X[,i-1] + (theta1 - theta2*X[,i-1])*dt + theta3*sqrt(X[,i-1])*Z[,i-1]
  }
  return(OU = list(X = X,
                    Z = Z))
}

Assetpaths = function(S0, K, sigma, T.t, r, nsim, nstep, V, Z.V){
  dt = T.t/nstep
  S = matrix(0, nsim, (1+nstep)) 
  S[,1] = rep(S0,nsim)
  for (i in 2:(nstep+1)){ 
    eps = rnorm(nsim)
    S[,i] = S[,i-1]*exp((r-0.5*V[,i-1])*dt + sqrt(V[,i-1])*(rho*Z.V[,i-1] + sqrt((1-rho^2))*dt*eps))
  }
  return(S)
}

# fast algo
expBes = function(x,nu){
  mu = 4*nu ^2
  A1 = 1
  A2 = A1 * (mu - 1) / (1 * (8*x))
  A3 = A2 * (mu - 9) / (2 * (8*x))
  A4 = A3 * (mu - 25) / (3 * (8*x))
  A5 = A4 * (mu - 49) / (4 * (8*x))
  A6 = A5 * (mu - 81) / (5 * (8*x))
  A7 = A6 * (mu -121) / (6 * (8*x))
  1/ sqrt(2*pi*x) * (A1 - A2 + A3 - A4 + A5 - A6 + A7)
}

dcCIR = function (x, t, x0 , theta , log = FALSE ){
  c = 2* theta[2] /((1 - exp(- theta[2] *t))* theta[3]^2)
  ncp = 2*c*x0*exp(- theta[2] *t)     # non centrality param of chi sq
  df = 4* theta[1] / theta[3]^2       # df of condition prob of the chi sq con. prob
  u = c*x0* exp (- theta [2] *t)
  v = c*x
  q = 2* theta [1] / theta [3]^2 -1
  lik = ( log (c) - (u+v) + q/2 * log (v/u) + log ( expBes ( 2* sqrt (u*v), q))
           + 2* sqrt (u*v))
  if(!log )
    lik = exp(lik)
  lik
}

CIR.lik = function(theta1 , theta2 , theta3 ) {
  n = length(X)
  dt = deltat(X)
  -sum(dcCIR(x=X[2: n], t=dt , x0=X[1:(n -1)] ,theta=c(theta1 , theta2 , theta3 ),
             log = TRUE ))
}

garchpara = garch(ret, order = c(1,1))

#### stochastic variance ####
## simulate path
d = 20                                                         # Trading days to maturity 
T.t = d/252                                                    # Annualized time-to-maturity
nstep = d*2                                                    # num of steps
nsim = 10^6                                                    # num of paths
theta2 = as.numeric(unlist(garchpara$coef[2]))                 # k
theta1 = as.numeric(unlist(garchpara$coef[1]))*theta2          # k*long_run_mean
theta3 = sqrt(2*theta1)-as.numeric(unlist(garchpara$coef[3]))  # volatility of the stock price volatility (2*k*theta > sigma^2)
x0 = 0.2                                                       # stock price volatility starting point

OU = ornstein_uhlenbeck(T.t, nstep, nsim, theta1, theta2, theta3, x0) 
V = OU$X #sim vola matrix

# diagnostic plot
par(mfrow=c(1,2))
matplot(t(V[1:1000,]), type='l', main = "Simulated Volatility Matrix", ylab = "Values", xlab = "Steps")
hist(V[,ncol(V)], breaks = 100, col='2', main = "Volatility Outcome Histogram", xlab = "Volatility Outcome")

## parameters calibration
X = V[1,] #Set fit Garch(1,1)
fit = mle(CIR.lik, start = list(theta1 = theta1, theta2 = theta2, theta3 = theta3),
           method ="BFGS")
fit

long_run_mean = fit@coef[1]/fit@coef[2]
k = fit@coef[2]
sigma = fit@coef[3]
2*k*long_run_mean > sigma^2 # TRUE

OU = ornstein_uhlenbeck(T.t, nstep, nsim, theta1 <- fit@coef[1], theta2 <- fit@coef[2], theta3 <- fit@coef[3], x0)
V = OU$X
Z.V = OU$Z

# diagnostic
par(mfrow=c(1,2),oma=c(0,0,2,0))
matplot(t(V[1:10000,]), type='l', main = "", ylab = "Values", xlab = "Steps")
hist(V[,ncol(V)], breaks = 100, col='blue', main = "", xlab = "Volatility Outcome")
title("Simulated Volatility Matrix with Calibration", outer=TRUE)


z = (V[,nstep+1] - mean(V[,nstep+1]) )/sd(V[,nstep+1])
par(mfrow=c(1,2))
hist(z, breaks = 100, main = "", xlab = "Outcomes")
hist(rnorm(nsim, 0, 1), breaks = 100, main="", xlab = "Outcomes")
title("Normalized Simulated Variance vs Variance of a Normal", outer = TRUE)

par(mfrow=c(1,1))
hist(z, breaks = 100,  prob=T)
lines(density(rnorm(nsim, 0, 1) ), col='2')

skewness(z) ; skewness(rnorm(nsim, 0, 1))
kurtosis(z) ; kurtosis(rnorm(nsim, 0, 1))

#### log(S_t) process ####
d = d                             # Trading days to maturity 
T.t = T.t                         # Annualized time-to-maturity
nstep = nstep                     # num of steps
nsim = nsim                       # num of paths
S0 = S0                           # S(t=0)
rho = 0.35                        # correlation(S,V)

assetpath = Assetpaths(S0, K, sigma, T.t, r, nsim, nstep, V, Z.V)

# diagnostic
par(mfrow=c(1,2))
matplot(t(assetpath[1:10000,]), type='l')
hist(assetpath[,ncol(assetpath)], breaks = 100, col='2')

z = (assetpath[,nstep+1] - mean(assetpath[,nstep+1]) )/sd(assetpath[,ncol(assetpath)])
hist(z, breaks = 100)
hist(rnorm(nsim, 0, 1), breaks = 100)

par(mfrow=c(1,1))
hist(z, breaks = 100,  prob=T)
lines(density(rnorm(nsim, 0, 1) ), col='2')
skewness(z) ; skewness(rnorm(nsim, 0, 1))
kurtosis(z) ; kurtosis(rnorm(nsim, 0, 1))


#### HESTON ####
S0 = S0                                              #Initial stock price
n = 10^6                                             #Number of simulations (feel free to reduce this)
freq = "monthly"                                     #Sampling frequency
kappa = unlist(as.numeric(fit@coef[2]))              #Volatility mean-reversion speed
volvol = unlist(as.numeric(fit@coef[3]))             #Volatility of volatility
rho = -0.5                                           #Correlation between stoch. vol and spot prices
V0 = 0.2                                             #Initial variance
theta = unlist(as.numeric(fit@coef[1]/fit@coef[2]))  #long-term variance
r0 = r                                               #Initial short rate

horizon <- T                                          # Options maturities
strikes <- K                                          # Options' exercise prices

# Simulation of shocks with given correlation
set.seed(1234) # reproducibility seed
shocks = simshocks(n =  n,
                    horizon =  horizon,
                    frequency =  freq,
                    method = "anti",
                    family = 1, par =  rho)

#  Stochastic volatility  simulation
sim.vol = simdiff(n =  n, horizon =  horizon,
                   frequency =  freq, model = "CIR", x0 =  V0,
                   theta1 =  kappa*theta, theta2 =  kappa,
                   theta3 =  volvol, eps =  shocks[[1]])

# Stock prices simulation
sim.price = simdiff(n = n, horizon = horizon,
                     frequency = freq, model = "GBM", x0 = S0,
                     theta1 = r0, theta2 = sqrt(sim.vol),
                     eps = shocks[[2]])


# Stock price at maturity
S_T = sim.price[nrow(sim.price), ]

### Monte Carlo prices
#### Estimated Monte Carlo price
discounted.payoff = function(x)
{
  (S_T - x)*(S_T - x > 0)*exp(-r0*horizon)
}
mcprices = sapply(strikes, function(x)
  mean(discounted.payoff(x)))

#### 95% Confidence interval around the estimation
mcprices95 = sapply(strikes,  function(x)
  t.test(discounted.payoff(x),
         conf.level = 0.95)$conf.int)

#### 'Analytical' prices given by 'callHestoncf'
pricesAnalytic = sapply(strikes, function(x)
  callHestoncf(S = S0, X = x, tau = horizon,
               r = r0, q = 0, v0 = V0, vT = theta,
               rho = rho, k = kappa, sigma = volvol))

HestonCall = data.frame(cbind(strikes, mcprices,
                            t(mcprices95), pricesAnalytic))
colnames(HestonCall) = c("strikes", "mcprices", "lower95",
                       "upper95", "pricesAnalytic")

print(HestonCall)




####  JUMP  ####
library(NMOF)
callMerton(S, K, T, r, 0, sigma, lambda, muJ, vJ, N, TRUE)

#Function disclosed
# S current stock price
# X strike price
# tau time to maturity
# r risk-free rate
# q dividend rate
# sigma variance
# lambda jump intensity
# muJ mean jump-size
# vJ variance of log jump-size
# N The number of jumps. See Details.
# implVol compute equivalent Black???Scholes???Merton volatility? Default is FALSE.

#Data
S = S0     
X = K
tau = T
r = r 
q = 0
sigma = sigma
lambda = 1;
N = 1000

#Count jumps and derive measures
####  JUMP COUNTER IN A TIME SERIES ####
#Data
jcount = NULL
jcount[1] = 0
for (i in 2:length(stockprice)){
  if (stockprice[i] > stockprice[i-1]*1.025 & stockprice[i] < stockprice[i-1]*1.05){
    jcount[i] = 1
  } else if (stockprice[i] >= stockprice[i-1]*1.5 & stockprice[i] < stockprice[i-1]*1.075 ){
    jcount[i] = 2
  } else if (stockprice[i] >= stockprice[i-1]*1.075 & stockprice[i] < stockprice[i-1]*1.10 ){
    jcount[i] = 3
  } else if (stockprice[i] >= stockprice[i-1]*1.10 & stockprice[i] < stockprice[i-1]*1.125 ){
    jcount[i] = 4
  } else if (stockprice[i] <= stockprice[i-1]*0.975 & stockprice[i] > stockprice[i-1]*0.95){
    jcount[i] = -1
  } else if (stockprice[i] <= stockprice[i-1]*0.95 & stockprice[i] > stockprice[i-1]*0.925){
    jcount[i] = -2
  } else if (stockprice[i] <= stockprice[i-1]*0.925 & stockprice[i] > stockprice[i-1]*0.875){
    jcount[i] = -3
  } else if (stockprice[i] <= stockprice[i-1]*0.875 & stockprice[i] > stockprice[i-1]*0.85){
    jcount[i] = -4
  } else if (stockprice[i] >= stockprice[i-1]*1.125 & stockprice[i] < stockprice[i-1]*1.15){
    jcount[i] = 5
  } else if (stockprice[i] >= stockprice[i-1]*1.15 & stockprice[i] < stockprice[i-1]*1.175 ){
    jcount[i] = 6
  } else if (stockprice[i] <= stockprice[i-1]*0.85 & stockprice[i] > stockprice[i-1]*0.825){
    jcount[i] = -5
  } else if (stockprice[i] <= stockprice[i-1]*0.825 & stockprice[i] > stockprice[i-1]*0.8){
    jcount[i] = -6
  } else {
    jcount[i] = 0
  }
}

plot(sort(jcount), main = "Jump Counter Outcome", ylab = "Jump Size", xlab= "Time")
muJ = mean(jcount)
logjcount = pmax(log(jcount), 0)
logjcount[is.na(logjcount)] = 0
vJ = var(logjcount)

#Pricing Formula
callBSM = function(S, X, tau, r, q, sigma) {
  d1 = (log(S/X) + (r - q + sigma/2) * tau)/(sqrt(sigma) * sqrt(tau))
  d2 = d1 - sqrt(sigma) * sqrt(tau)
  S * exp(-q * tau) * pnorm(d1) - X * exp(-r * tau) * pnorm(d2)
}

lambda2 = lambda * (1 + muJ)

JumpCall = 0
for (n in 0:N) {
  v_n = sigma + n * vJ/tau
  r_n = r - lambda * muJ + n * log(1 + muJ)/tau
  JumpCall = JumpCall + (exp(-lambda2 * tau) * (lambda2 * tau)^n) * callBSM(S, X, tau, r_n, q, v_n)/factorial(n)
}
JumpCall



####  FRACTAL ####
#R/S ANALYSIS
#### HURST PRE-OPT PRICING ####

library(tseries)

#Data
prices = stockprice
#prices
plot(prices, type='l')
head(prices)

#log-returns
return = diff(log(prices))
summary(return)
plot(return, type = 'l')

#splits
split(return, c(1:512), drop = FALSE)

#decomposition loop
k = 9
for(i in 1:k){
  h = i-1
  j = 2^h
  assign(paste("split",i,sep=""),(split(return, ceiling(seq_along(return)/(length(return)/j)))))
}

split_one = list(split1,split2,split3,split4,split5,split6,split7,split8,split9)
avg.vec = NULL
stdev.vec = NULL
for(i in 1:length(split_one)){
  avg = NULL
  stdev = NULL
  for(j in 1:length(split_one[[i]])){
    a = mean(unlist(split_one[[i]][j]))
    avg[j] = a 
    b = sd(unlist(split_one[[i]][j])) 
    stdev[j] = b
  }
  avg.vec = append(avg.vec, avg)
  stdev.vec = append(stdev.vec, stdev)
}




#errors' list
##########################################################
lappend = function(list, object) {
  list[[length(list)+1]] = object
  return(list)
}
##########################################################


list.err = list()
counter = 0

for(i in 1:length(split_one)){
  errors = list()
  for(j in 1:length(split_one[[i]])){
    
    a = unlist(split_one[[i]][j])-avg.vec[j+counter]
    errors[j] = list(cumsum(a))
    
  }
  list.err[[i]] = errors
  counter <- j + counter
  
}
names(vec) = NULL



#min & max error
min.vec = NULL
max.vec = NULL
for(i in 1:length(list.err)){
  min = NULL
  max = NULL
  for(j in 1:length(list.err[[i]])){
    c = min(unlist(list.err[[i]][j]))
    min[j] = c 
    d = max(unlist(list.err[[i]][j]))
    max[j] = d
  }
  min.vec = append(min.vec,min)
  max.vec = append(max.vec,max)
}

#Range vector
range = max.vec - min.vec

#Rescaled range
rescaled = range/stdev


###############################################################
i = 0
lis = list()
for(i in 1:9){
  lis = split(avg.vec, sort(rank(avg.vec)/ 2^i))
}
###############################################################

# group rescaled range (R/S)

resc = list()
j = 0
split = 8
for( i in 0:split) {
  k = 2^i
  j = k + j
  resc[[i+1]] = rescaled[c(k:j)]
}


#R/S groups avg

avg.vec.resc = NULL
avg.resc = NULL
for(i in 1:length(resc)){
  h = mean(unlist(resc[[i]]))
  avg.resc[i] = h 
  
  avg.vec.resc[i] = rbind(avg.resc[i])
}

avg.vec.resc

#N
n = NULL
for(i in 1:length(split_one)){
  number = NULL
  for(j in 1:length(split_one[[i]])){
    g = length(unlist(split_one[[i]][j]))
    number[j] = g
  }
  n<-append(n,number)
}


n_0 = list()
g = 0
split = 8
for( y in 0:split) {
  u = 2^y
  g = u + g
  n_0[[y+1]] = n[c(u:g)]
}
n.vec = NULL
n.avg = NULL
for(i in 1:length(n_0)){
  h = mean(unlist(n_0[[i]]))
  n.avg[i] = h 
  
  n.vec[i] = rbind(n.avg[i])
}

n.vec


obs = log(n.vec)
R_S = log(avg.vec.resc)

#OLS
reg = (lm(R_S~obs))
summary(reg)

H = as.numeric(reg$coefficients[2])

#PRICING
FBMCall = function(S0,K,sigma,tau,r,H){
  dhat1 = (log(S0/K) + r*tau + 0.5*sigma^2*tau^(2*H))/(sigma*tau^H)
  dhat2 = dhat1 - sigma*tau^H
  FBMCall = S0 * pnorm(dhat1) - K*exp(-r*tau) * pnorm(dhat2)
  return(FBMCall)
}
A = FBMCall(S0,K,sigma,tau,r,H)


BScall=function(S0,K,sigma,T,r){
  d1=(log(S0/K)+(r+0.5*sigma^2)*T)/(sigma*sqrt(T))
  d2=d1-sigma*sqrt(T);
  Price=S0*pnorm(d1)-K*exp(-r*T)*pnorm(d2)
  return(Price)
}

B = BScall(S0,K,sigma,T,r)

B-A
#plot(B-A)


#### CHECK ####

Call
Cbsmc
Cbs
Cbsmcav
Cbsmccv
#######
HestonCall
#######
JumpCall
#######
A

OptionPricesAAPL = data.frame(Call@price, Cbs, Cbsmc[[1]], Cbsmcav[[1]], Cbsmccv[[1]], HestonCall$pricesAnalytic, JumpCall, A)
colnames(OptionPricesAAPL) = c("BSbi", "BSf", "BSMC", "BSMCav", "BSMCcv", "HESTON", "JUMP", "FBS")
rownames(OptionPricesAAPL) = "Call Price"
View(OptionPricesAAPL)


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
