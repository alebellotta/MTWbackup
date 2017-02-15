library("somebm")
FBM1 = fbm(hurst = 0.1, n = 10000)
plot(FBM1)

FBM5 = fbm(hurst = 0.5, n = 10000)
plot(FBM5)

FBM10 = fbm(hurst = 0.9, n = 10000)
plot(FBM10)

par(mfrow=c(1,3))
plot(FBM1, main = "H = 0.1 - Negative Correlation", ylab = "Value")
plot(FBM5, main = "H = 0.5 - Brownian Motion", ylab = "Value")
plot(FBM10, main = "H = 0.9 - Positive Correlation", ylab = "Value")
