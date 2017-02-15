####  DATA  ####
#get risk free rate#
getSymbols.FRED("TB1YR", env = globalenv())
TB1YR = as.data.frame(TB1YR)
TB1YR = TB1YR$TB1YR
r = tail(TB1YR, n=1)

#get data#
ticker = "AAPL"
stockprice = as.numeric(get.hist.quote(instrument = ticker, quote = "Close"))

ret = diff(log(stockprice))
S0 = tail(stockprice, n=1) 


K = S0*1.2  #c(S0*0.2, S0*0.4, S0*0.6, S0*0.8, S0, S0*1.2, S0*1.4, S0*1.6, S0*1.8, S0*2)
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
return<- diff(log(prices))
summary(return)
plot(return, type = 'l')

#splits
split(return, c(1:512), drop = FALSE)

#decomposition loop
k <- 9
for(i in 1:k){
  h <- i-1
  j <- 2^h
  assign(paste("split",i,sep=""),(split(return, ceiling(seq_along(return)/(length(return)/j)))))
}

split_one<- list(split1,split2,split3,split4,split5,split6,split7,split8,split9)
avg.vec<- NULL
stdev.vec<- NULL
for(i in 1:length(split_one)){
  avg <- NULL
  stdev<-NULL
  for(j in 1:length(split_one[[i]])){
    a<-mean(unlist(split_one[[i]][j]))
    avg[j]<-a 
    b<-sd(unlist(split_one[[i]][j])) 
    stdev[j]<-b
  }
  avg.vec<- append(avg.vec, avg)
  stdev.vec<- append(stdev.vec, stdev)
}




#errors' list
##########################################################
lappend <- function(list, object) {
  list[[length(list)+1]] <- object
  return(list)
}
##########################################################


list.err<- list()
counter <- 0

for(i in 1:length(split_one)){
  errors <- list()
  for(j in 1:length(split_one[[i]])){
    
    a <- unlist(split_one[[i]][j])-avg.vec[j+counter]
    errors[j] <- list(cumsum(a))
    
  }
  list.err[[i]] <- errors
  counter <- j + counter
  
}
names(vec)<- NULL



#min & max error
min.vec <- NULL
max.vec<- NULL
for(i in 1:length(list.err)){
  min <- NULL
  max <- NULL
  for(j in 1:length(list.err[[i]])){
    c<-min(unlist(list.err[[i]][j]))
    min[j]<-c 
    d<- max(unlist(list.err[[i]][j]))
    max[j]<-d
  }
  min.vec <- append(min.vec,min)
  max.vec <- append(max.vec,max)
}

#Range vector
range<- max.vec - min.vec

#Rescaled range
rescaled <- range/stdev


###############################################################
i<-0
lis<-list()
for(i in 1:9){
  lis<-split(avg.vec, sort(rank(avg.vec)/ 2^i))
}
###############################################################

# group rescaled range (R/S)

resc <- list()
j <- 0
split <- 8
for( i in 0:split) {
  k <- 2^i
  j <- k + j
  resc[[i+1]] <- rescaled[c(k:j)]
}


#R/S groups avg

avg.vec.resc<- NULL
avg.resc <- NULL
for(i in 1:length(resc)){
  h<-mean(unlist(resc[[i]]))
  avg.resc[i]<-h 
  
  avg.vec.resc[i]<- rbind(avg.resc[i])
}

avg.vec.resc

#N
n<-NULL
for(i in 1:length(split_one)){
  number<-NULL
  for(j in 1:length(split_one[[i]])){
    g <- length(unlist(split_one[[i]][j]))
    number[j] <- g
  }
  n<-append(n,number)
}


n_0<- list()
g <- 0
split <- 8
for( y in 0:split) {
  u <- 2^y
  g <- u + g
  n_0[[y+1]] <- n[c(u:g)]
}
n.vec<- NULL
n.avg <- NULL
for(i in 1:length(n_0)){
  h<-mean(unlist(n_0[[i]]))
  n.avg[i]<-h 
  
  n.vec[i]<- rbind(n.avg[i])
}

n.vec


obs<-log(n.vec)
R_S<-log(avg.vec.resc)

#OLS
reg<-(lm(R_S~obs))
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
S <- S0     
X <- K
tau <- T
r <- r 
q <- 0
sigma = sigma
lambda <- 1;
N <- 75

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

plot(sort(jcount))
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
  JumpCall = JumpCall + (exp(-lambda2 * tau) * (lambda2 * tau)^n) * FBMCall(S, X, tau, r_n, q, v_n)/factorial(n)
}
JumpCall

