module MathToolsMod
using LeafPhotosynthesisMod:fluxes
using Leaf:leaf_params

export quadratic,beta_function,log_gamma_function, hybrid, zbrent, ∫, Avogadro, e2phot, ephoton
# hybrid             ! Solve for the root of a function using secant and Brent's methods
# zbrent             ! Use Brent's method to find the root of a function
# quadratic          ! Solve a quadratic equation for its two roots
# tridiag            ! Solve a tridiagonal system of equations
# beta_function      ! Evaluate the beta function at p and q: B(p,q)
# log_gamma_function ! Evaluate the log natural of the gamma function at x: ln(G(x))

const Avogadro = 6.022140857e23
const h = 6.62607004e-34          # [J s]         Planck's constant
const c = 299792458.0             #  [m s-1]       speed of light
#-----------------------------------------------------------------------
function quadratic(a, b, c)
  discr = b^2 - 4*a*c
  discr >= 0 ?   ( (-b + sqrt(discr))/(2a), (-b - sqrt(discr))/(2a) ) : 0.0 #error("Only complex roots")
end # function

function e2phot(lambda,E)
    #molphotons = e2phot(lambda,E) calculates the number of moles of photons
    #corresponding to E Joules of energy of wavelength lambda (m)
    e           = ephoton(lambda);
    photons     = E./e;
    return photons./Avogadro;
end

function ephoton(lambda)
    #E = phot2e(lambda) calculates the energy content (J) of 1 photon of
    #wavelength lambda (m)
    return (h*c)./lambda;           # [J]           energy of 1 photon
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

function hybrid(flux::fluxes,leaf::leaf_params, func::Function, xa::Number, xb::Number, tol::Number)
    #
    # !DESCRIPTION:
    # Solve for the root of a function, given initial estimates xa and xb.
    # The root is updated until its accuracy is tol.
    #
    # !USES:
    #
    # !ARGUMENTS:
    #procedure (xfunc) :: func         ! Function to solve
    #real(r8), intent(in) :: xa, xb    ! Initial estimates of root
    #real(r8), intent(in) :: tol       ! Error tolerance
    #type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    itmax = 40
    x0 = x1 = f0 = f1 = 0.0
    x0 = xa
    f0 = func(x0,flux, leaf)
    if (f0 == 0.0) return  x0 end

    x1 = xb
    f1 = func(x1,flux, leaf)
    if (f1 == 0.0) return x1 end

    if (f1 < f0)
       minx = x1
       minf = f1
    else
       minx = x0
       minf = f0
    end

    # First use the secant method, and then use the brent method as a backup
    for iter = 0:itmax
       iter = iter + 1
       dx = -f1 * (x1 - x0) / (f1 - f0)
       x = x1 + dx
       if (abs(dx) < tol)
          x0 = x
          break
       end
       x0 = x1
       f0 = f1
       x1 = x
       f1 = func(x1,flux, leaf)
       if (f1 < minf)
          minx = x1
          minf = f1
       end

       # If a root zone is found, use the brent method for a robust backup strategy
       if (f1 * f0 < 0.0)
          x = zbrent(flux,  leaf,func,x0, x1, xtol=tol)
          x0 = x
          break
       end

       # In case of failing to converge within itmax iterations stop at the minimum function
       if (iter > itmax)
          f1 = func(minx,flux, leaf)
          x0 = minx
          break
       end

    end

    return x0
end #function hybrid

function zbrent(flux::fluxes,leaf::leaf_params,f::Function, x0::Number, x1::Number, args::Tuple=();
               xtol::AbstractFloat=1e-7, ytol=2eps(Float64),
               maxiter::Integer=50)
    EPS = eps(Float64)
    y0 = f(x0,flux,leaf)
    y1 = f(x1,flux,leaf)
    if abs(y0) < abs(y1)
        # Swap lower and upper bounds.
        x0, x1 = x1, x0
        y0, y1 = y1, y0
    end
    x2 = x0
    y2 = y0
    x3 = x2
    bisection = true
    for _ in 1:maxiter
        # x-tolerance.
        if abs(x1-x0) < xtol
            return x1
        end

        # Use inverse quadratic interpolation if f(x0)!=f(x1)!=f(x2)
        # and linear interpolation (secant method) otherwise.
        if abs(y0-y2) > ytol && abs(y1-y2) > ytol
            x = x0*y1*y2/((y0-y1)*(y0-y2)) +
                x1*y0*y2/((y1-y0)*(y1-y2)) +
                x2*y0*y1/((y2-y0)*(y2-y1))
        else
            x = x1 - y1 * (x1-x0)/(y1-y0)
        end

        # Use bisection method if satisfies the conditions.
        delta = abs(2EPS*abs(x1))
        min1 = abs(x-x1)
        min2 = abs(x1-x2)
        min3 = abs(x2-x3)
        if (x < (3x0+x1)/4 && x > x1) ||
           (bisection && min1 >= min2/2) ||
           (!bisection && min1 >= min3/2) ||
           (bisection && min2 < delta) ||
           (!bisection && min3 < delta)
            x = (x0+x1)/2
            bisection = true
        else
            bisection = false
        end

        y = f(x,flux,leaf)
        # y-tolerance.
        if abs(y) < ytol
            return x
        end
        x3 = x2
        x2 = x1
        if sign(y0) != sign(y)
            x1 = x
            y1 = y
        else
            x0 = x
            y0 = y
        end
        if abs(y0) < abs(y1)
            # Swap lower and upper bounds.
            x0, x1 = x1, x0
            y0, y1 = y1, y0
        end
    end
    error("Max iteration exceeded")
end

function ∫(x::Vector, f::Vector)
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

end # module MathToolsMod
