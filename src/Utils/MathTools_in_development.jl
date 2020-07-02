###############################################################################
#
# Calculate the frequency of leaf inclination in a canopy
# These functions passed the FT tests
# Documentation is not yet finished
#
###############################################################################
#=
"""
    dcum(a::FT, b::FT, t::FT)

TODO Add function description
"""
=#
function dcum(a::FT, b::FT, t::FT) where {FT}
    y = FT(0.0)
    if a >= 1
        f = 1 - cosd(t)
    else
        epsi = 1e-8
        delx = FT(1.0)
        x    = 2 * deg2rad(t)
        p    = x
        # add a iteration number (n) control to avoid rounding issue in Float32
        n    = 0
        while (delx >= epsi) && (n<50)
            n   += 1
            y    = a*sin(x) + b*sin(2x)/2
            dx   = (y - x + p) / 2
            x   += dx
            delx = abs(dx)
        end
    	f = (2*y + p) / pi
    end
    return f
end



#=
"""
    dladgen(a::FT, b::FT, litab_bnd::Array)

TODO Calculate the freqency of WHAT?
"""
=#
function dladgen(a::FT, b::FT, litab_bnd::Array{FT,2}) where {FT}
    @assert a+b< 1
    freq = similar(litab_bnd[:,1])
    for i in eachindex(freq)
        freq[i] = dcum(a, b, litab_bnd[i,2]) - dcum(a, b, litab_bnd[i,1])
    end
    return freq
end








###############################################################################
#
# psofunction
# This function passed the FT test
# Documentation is not yet finished
#
###############################################################################
#=
"""
    psofunction(K::FT, k::FT, Ω::FT, LAI::FT, q::FT, dso::FT, xl::FT)

TODO explain the variables
Return the probability of observing a sunlit leaf at depth `xl` (`pso`, see eq 31 in vdT 2009), given
- `xl` Leaf depth in the canopy
"""
=#
function psofunction(K::FT, k::FT, Ω::FT, LAI::FT, q::FT, dso::FT, xl::FT) where {FT}
    # [nl+1]  factor for correlation of Ps and Po
    if dso != 0
        alf = dso / q * 2 / (k + K)
        return exp( (K+k)*Ω*LAI*xl + sqrt(K*k)*Ω*LAI / alf * (1 - exp(xl*alf)) )
    else
        return exp( (K + k -sqrt(K*k)) * Ω * LAI * xl )
    end
end










function lower_quadratic2!(a, b, c,var)
    #r = similar(a)
    discr = b.^2 .- 4 * a.*c
    #@show discr
    var[discr .>= 0] = (-b .- sqrt.(discr[discr .>= 0]))./2a
    #a[discr .< 0].=NaN
  end # function

  function lower_quadratic!(a, b, c,vari)
    discr = b^2 - 4*a*c
    if discr >= 0
        vari=(-b - sqrt(discr))/2a;
    else
        vari=NaN
    end
  end # function

  function ephoton(λ::Array)
      FT = eltype(λ)
      #E = phot2e(lambda) calculates the energy content (J) of 1 photon of
      #wavelength lambda (m)
      return (FT(H_PLANCK)*FT(LIGHT_SPEED))./λ;           # [J]           energy of 1 photon
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


  #=
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
  =#
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
