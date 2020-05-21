

# hybrid             ! Solve for the root of a function using secant and Brent's methods
# zbrent             ! Use Brent's method to find the root of a function
# quadratic          ! Solve a quadratic equation for its two roots
# tridiag            ! Solve a tridiagonal system of equations
# beta_function      ! Evaluate the beta function at p and q: B(p,q)
# log_gamma_function ! Evaluate the log natural of the gamma function at x: ln(G(x))

using IncGammaBeta


export quadratic, lower_quadratic,
       beta_function,
       log_gamma_function,
       ∫,
       fast∫,
       Avogadro,
       e2phot,
       ephoton,
       Weibull,
       IntWeibull

# hybrid             ! Solve for the root of a function using secant and Brent's methods
# zbrent             ! Use Brent's method to find the root of a function
# quadratic          ! Solve a quadratic equation for its two roots
# tridiag            ! Solve a tridiagonal system of equations
# beta_function      ! Evaluate the beta function at p and q: B(p,q)
# log_gamma_function ! Evaluate the log natural of the gamma function at x: ln(G(x))

const Avogadro = 6.022140857e23
const h = 6.62607004e-34          # [J s]         Planck's constant
const c = 299792458.0             #  [m s-1]       speed of light

function quadratic(a, b, c)
  discr = b^2 - 4*a*c
  discr >= 0 ?   ( (-b + sqrt(discr))/2a, (-b - sqrt(discr))/2a ) : error("imaginary roots in quadratic")
end # function

function lower_quadratic(a, b, c)
  discr = b^2 - 4*a*c
  discr >= 0 ?   (-b - sqrt(discr))/2a  : error("imaginary roots in quadratic")
end # function

function e2phot(λ::Array,E::Array)
    FT = eltype(E)
    #molphotons = e2phot(lambda,E) calculates the number of moles of photons
    #corresponding to E Joules of energy of wavelength lambda (m)
    fac = FT(h)*FT(c)*FT(Avogadro)
    (E.*λ)*(FT(1e-9)/fac)
end

function ephoton(λ::Array)
    FT = eltype(λ)
    #E = phot2e(lambda) calculates the energy content (J) of 1 photon of
    #wavelength lambda (m)
    return (FT(h)*FT(c))./λ;           # [J]           energy of 1 photon
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


"""
Simpson rule for irregularly spaced data.

    Parameters
    ----------
    x : list or np.array of floats
    Sampling points for the function values
    f : list or np.array of floats
    Function values at the sampling points

    Returns
    -------
    float : approximation for the integral
"""
function ∫(x::Vector, f::Vector)
    N = length(x)
    h = diff(x)
    #println(length(x)," ", length(f))
    result = typeof(x[1])(0.0)
    for i=2:2:N-1
    hph = h[i] + h[i - 1]
    result += f[i] * ( h[i]^3 + h[i - 1]^3 + 3. * h[i] * h[i - 1] * hph )/ ( 6 * h[i] * h[i - 1] )
    result += f[i - 1] * ( 2 * h[i - 1]^3 - h[i]^3+ 3. * h[i] * h[i - 1]^2) / ( 6 * h[i - 1] * hph)
    result += f[i + 1] * ( 2 * h[i]^3 - h[i - 1]^3 + 3. * h[i - 1] * h[i]^2)/ ( 6 * h[i] * hph )
    end

    if N % 2 == 0
    result += f[N] * ( 2 * h[N - 1]^2+ 3 * h[N - 2] * h[N - 1]) / ( 6 * ( h[N - 2] + h[N - 1] ) )
    result += f[N - 1] * ( h[N - 1]^2+ 3*h[N - 1]* h[N - 2] ) / ( 6 * h[N - 2] )
    result -= f[N - 2] * h[N - 1]^3/ ( 6 * h[N - 2] * ( h[N - 2] + h[N - 1] ) )
    end
    return result
end

# Already uses diff here
function fast∫(dx::Vector, f::Vector)
    N = length(f)
    result = eltype(f)(0.0)
    @inbounds for i=1:N
      result += f[i]* dx[i]
    end
    return result
end


#-----------------------------------------------------------------------
# integral of weibull function for conductances
function IntWeibull(psis,psil,psil50,ck) # set hydraulic conductivity
# checked - this is correct
    one_c = 1.0/ck;
    b1 = log(2)*(psil/psil50)^ck;
    b2 = log(2)*(psis/psil50)^ck;
    integral = psil50/(ck*log(2)^one_c) * ( inc_gamma_upper(one_c,b1) - inc_gamma_upper(one_c,b2) )
    #println("Int weibull= ",integral)
    return max(integral,0.0);

    # -(psi_50*gamma_incomplete(1/c,(ln(2)*x^c)/psi_50^c))/(c*ln(2)^(1/c))
end

# integral of weibull function for conductances
function Weibull(psil,psil50,ck) # set hydraulic conductivity
    return 2^(-(psil/psil50)^ck) ;
end
