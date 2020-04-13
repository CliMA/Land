####
#### Math tools
####

#-----------------------------------------------------------------------
# !DESCRIPTION:
# Math tools
#
# !USES:
#use shr_kind_mod, only : r8 => shr_kind_r8
#use clm_varctl, only : iulog
#use abortutils, only : endrun
#use CanopyFluxesMultilayerType, only : mlcanopy_type
#!
#! !PUBLIC TYPES:
#implicit none
#!
#! !PUBLIC MEMBER FUNCTIONS:
export quadratic,beta_function,IntWeibull,log_gamma_function
# hybrid             ! Solve for the root of a function using secant and Brent's methods
# zbrent             ! Use Brent's method to find the root of a function
# quadratic          ! Solve a quadratic equation for its two roots
# tridiag            ! Solve a tridiagonal system of equations
# beta_function      ! Evaluate the beta function at p and q: B(p,q)
# log_gamma_function ! Evaluate the log natural of the gamma function at x: ln(G(x))

export quadratic, IntWeibull, beta_function, log_gamma_function

#-----------------------------------------------------------------------
function quadratic(a, b, c)
  discr = b^2 - 4*a*c
  discr >= 0 ?   ( (-b + sqrt(discr))/(2a), (-b - sqrt(discr))/(2a) ) : error("Only complex roots")
end # function

#-----------------------------------------------------------------------
# integral of weibull function for conductances
function IntWeibull(psis,psil,psil50,ck) # set hydraulic conductivity
    return 2^(-(psis/psil50)^ck) - 2^(-(psil/psil50)^ck);
end

#-----------------------------------------------------------------------
function beta_function(p, q)
#
# !DESCRIPTION:
# Return the value of the beta function evaluated at p and q: B(p,q)
#
# !USES:
#
exp(log_gamma_function(p) + log_gamma_function(q) - log_gamma_function(p+q))

end #function beta_function

#-----------------------------------------------------------------------
function log_gamma_function(x)
#
# !DESCRIPTION:
# Return the value of the log natural of the gamma function evaluated at x: ln(G(x))
#
  coef = [76.18009172947146, -86.50532032941677,24.01409824083091, -1.231739572450155, 0.1208650973866179e-02, -0.5395239384953e-05];
  stp = 2.5066282746310005;
  y = x
  tmp = x + 5.5
  tmp = (x + 0.5) * log(tmp) - tmp
  ser = 1.000000000190015;
  for co in coef
    y = y + 1.0
    ser = ser + co / y
  end
  return tmp + log(stp*ser/x)
end # function log_gamma_function
