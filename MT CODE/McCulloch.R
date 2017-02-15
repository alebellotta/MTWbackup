#### McCulloch Option Pricing ####
stable = rstable(length(ret), a, b, c, d)
norm = rnorm(length(ret), 0, 1)
plot(sort(stable))
plot(sort(norm))

H = 0.540292

a = 1/H
b = 1
c = 1
d = 0

S = tail(AAPL, n=1)
X = 110
T = 1
r = r

c.one = (0.5*(1+b))^(1/a) * c
c.two = (0.5*(1-b))^(1/a) * c

FTR = S*exp(r*T)


#Compute integrals
#define integrands
integrand.one = function(x) {exp(-c.two*x^2) * stable[1] * (1-((c.two*x - log(FTR/X) + b*c^a * sec(0.5*pi*a))/c.one))}
integrand.two = function(x) {exp(-c.one*x^2) * stable[1] * ((c.two*x - log(FTR/X) + b*c^a * sec(0.5*pi*a))/c.one)}

#integrate
I.one = integrate(integrand.one, lower=0, upper=1, stop.on.error = FALSE)$value
I.two = integrate(integrand.two, lower = 0, upper = 1)$value

C = FTR * exp(-r*T + c.two^a*sec(0.5*pi*a)) * I.one
  - X * exp(-r*T + c.one^a*sec(0.5*pi*a)) * I.two

for (i in 1:length(stable)){
  int.one[i+1] = int.one[i]
}

int = function(x){exp}






# Sample data
x <- tail(ret, n=100)

# Density (I reparametrize it to remove the constraints 
# on the parameters)
library(fBasics)
f <- function(u,a,b,c,d) {
  cat(a,b,c,d,"\n")  # Some logging (it is very slow)
  dstable(u, 2*exp(a)/(1+exp(a)), 2*exp(b)/(1+exp(b))-1, exp(c), d)
}

# Fit the distribution
library(MASS)
r <- fitdistr(x, f, list(a=1, b=0, c=log(mad(x)), d=median(x)))
r

# Graphical check
plot(
  qstable(ppoints(100),
          2*exp(r$estimate[1])/(1+exp(r$estimate[1])),
          2*exp(r$estimate[2])/(1+exp(r$estimate[2]))-1,
          exp(r$estimate[3]),
          r$estimate[4]
  ),
  sort(x)
)
abline(0,1)

# 
# library(fda)
# rangeval <- c(-3,3)#  set up range for density
# x <- x #  set up some standard normal data
# x[x < -3] <- -2.99#  make sure values within the range
# x[x >  3] <-  2.99
# 
# basisobj <- create.bspline.basis(rangeval, 11)#  set up basis for W(x)
# Wfd0 <- fd(matrix(0,11,1), basisobj)#  set up initial value for Wfdobj
# 
# WfdParobj <- fdPar(Wfd0)
# denslist <- density.fd(x, WfdParobj)#  estimate density
# 
# xval <- seq(-3,3,.2)#  plot density
# wval <- eval.fd(xval, denslist$Wfdobj)
# pval <- exp(wval)/denslist$C
# plot(xval, pval, type="l", ylim=c(0,0.4))
# points(x,rep(0,length(x)))






#### TEST FROM PAPER ####
theta = (pi * a)/2
lambda = 1
k = exp(lambda*d + lambda^a * c^a * 1/cos(theta))

CF = exp(1i * d * x + c^a * 1/cos(theta) * (lambda^a - (lambda - 1i*x)^a))

integrand.one = function(x) {exp(-c.two*x^2) * CF * (1-((c.two*x - log(FTR/X) + b*c^a * (1/cos(0.5*pi*a))/c.one)))}
integrand.two = function(x) {exp(-c.one*x^2) * CF * ((c.two*x - log(FTR/X) + b*c^a * (1/cos(0.5*pi*a))/c.one)}

integrate(integrand.one, 1, -1)









