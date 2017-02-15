#### callMerton Function Disclosed ####
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
S <- 100; X <- 100; tau <- 1
r <- 0.0075; q <- 0.00
sigma = 0.2
lambda <- 1;
N <- 1000

#Draw from random Poisson Distribution
poi = rpois(N, lambda)
muJ = mean(poi)
logpoi = pmax(log(poi), 0)
vJ = var(logpoi)

callBSM <- function(S, X, tau, r, q, sigma) {
    d1 <- (log(S/X) + (r - q + sigma/2) * tau)/(sqrt(sigma) * sqrt(tau))
    d2 <- d1 - sqrt(sigma) * sqrt(tau)
    S * exp(-q * tau) * pnorm(d1) - X * exp(-r * tau) * 
      pnorm(d2)
}
      
lambda2 <- lambda * (1 + muJ)

result <- 0
for (n in 0:N) {
    v_n <- sigma + n * vJ/tau
    r_n <- r - lambda * muJ + n * log(1 + muJ)/tau
    result <- result + (exp(-lambda2 * tau) * (lambda2 * tau)^n) * callBSM(S, X, tau, r_n, q, v_n)/factorial(n)
  }
  
# if (implVol) {
#     diffPrice <- function(vol, call, S, X, tau, r, q) {
#       d1 <- (log(S/X) + (r - q + vol^2/2) * tau)/(vol * 
#                                                     sqrt(tau))
#       d2 <- d1 - vol * sqrt(tau)
#       callBSM <- S * exp(-q * tau) * pnorm(d1) - X * exp(-r * 
#                                                            tau) * pnorm(d2)
#       call - callBSM
#     }
#     impliedVol <- uniroot(diffPrice, interval = c(1e-04, 
#                                                   2), call = result, S = S, X = X, tau = tau, r = r, 
#                           q = q)[[1L]]
#     result <- list(value = result, impliedVol = impliedVol)
#   }
#   result
# }

result 
callBSM(S, X, tau, r, q, sigma)
