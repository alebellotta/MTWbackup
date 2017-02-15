# Copyright (c) 2010, Dale Roberts.
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#   
#   * Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
# 
# * Neither the name of the author nor the names of its contributors
# may be used to endorse or promote products derived from this
# software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Dale Roberts <dale.roberts@anu.edu.au>
##


HestonInversionFormula = function(lambda, vbar, eta, rho, v0, r, tau, S0, K) {
    PIntegrand = function(u, lambda, vbar, eta, rho, v0, r, tau, S0, K, j) {
      F <- S0*exp(r*tau)
      x <- log(F/K)
      a <- lambda * vbar
      
      if (j == 1) {
        b = lambda - rho* eta
        alpha = - u^2/2 - u/2 * 1i + 1i * u
        beta = lambda - rho * eta - rho * eta * 1i * u
      } else {
        b = lambda
        alpha = - u^2/2 - u/2 * 1i
        beta = lambda - rho * eta * 1i * u
      }
      
      gamma = eta^2/2
      d = sqrt(beta^2 - 4*alpha*gamma)
      rplus = (beta + d)/(2*gamma)
      rminus = (beta - d)/(2*gamma)
      g = rminus / rplus
      
      D = rminus * (1 - exp(-d*tau))/(1-g*exp(-d*tau))
      C = lambda * (rminus * tau - 2/(eta^2) * log( (1-g*exp(-d*tau))/(1-g) ) )
      
      numerator = exp(C*vbar + D*v0 + 1i*u*x)
      denominator = (1i * u)
      Re(numerator/denominator)
    }
    
    P = function(lambda, vbar, eta, rho, v0, r, tau, S0, K, j) {
      value = integrate(PIntegrand, lower = 0, upper = Inf,
                         lambda, vbar, eta, rho, v0, r, tau,
                         S0, K, j, subdivisions=1000)$value
      0.5 + 1/pi * value
    }
    
    A = S0*P(lambda, vbar, eta, rho, v0, r, tau, S0, K, 1)
    B = K*exp(-r*tau)*P(lambda, vbar, eta, rho, v0, r, tau, S0, K, 0)
    A-B
}


HestonInversionFormula(lambda, vbar, eta, rho, v0, r, tau, S0, K)
