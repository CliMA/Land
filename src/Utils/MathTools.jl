module MathTools

using IncGammaBeta
using Parameters

using ..LandParameters

@unpack AVOGADRO,
        H_PLANCK,
        LIGHT_SPEED = LandParameters

export dladgen,
       e2phot,
       fast∫,
       lower_quadratic,
       psofunction,
       quadratic,
       volscatt




###############################################################################
#
# Convert energy to photons
# This function passed the FT test
# This function is documented in the Utils page
#
###############################################################################
_FAC = 1e-9 / (H_PLANCK * LIGHT_SPEED * AVOGADRO)

"""
    e2phot(λ::Array,E::Array)

Calculates the number of moles of photons, given
- `λ` An array of wave length in `[nm]`, converted to `[m]` by _FAC
- `E` Joules of energy
"""
function e2phot(λ::Array, E::Array)
    FT = eltype(E)
    fac::FT = _FAC
    return (E .* λ) .* fac
end








###############################################################################
#
# Integral functions
# This function passed the FT test
# This function is documented in the Utils page
#
###############################################################################
"""
    fast∫(dx::Vector, f::Vector)

A fast way of integrating functions, given
- `dx` Delta x for each x
- `f` f(x) for each x
"""
function fast∫(dx::Vector, f::Vector)
    FT = eltype(dx)
    N  = length(f)
    result::FT = 0.0
    @inbounds for i=1:N
      result += f[i]* dx[i]
    end
    return result
end








###############################################################################
#
# Quadratic function solver
# These functions passed the FT test
# THese functions area documented in the Utils page
#
###############################################################################
"""
lower_quadratic

Return the lower quadratic solution or NaN, given
- `a` Parameter in `a*x^2 + b*x + c = 0`
- `b` Parameter in `a*x^2 + b*x + c = 0`
- `c` Parameter in `a*x^2 + b*x + c = 0`
"""
function lower_quadratic(a::FT, b::FT, c::FT) where {FT}
    discr = b^2 - 4*a*c
    discr >= 0 ?   (-b - sqrt(discr))/2a : NaN
end




"""
    quadratic(a::FT, b::FT, c::FT)

Return a list of quadratic function solutions or an error message, given
- `a` Parameter in `a*x^2 + b*x + c = 0`
- `b` Parameter in `a*x^2 + b*x + c = 0`
- `c` Parameter in `a*x^2 + b*x + c = 0`
"""
function quadratic(a::FT, b::FT, c::FT) where {FT}
    discr = b^2 - 4*a*c
    discr >= 0 ? ( (-b + sqrt(discr))/2a, (-b - sqrt(discr))/2a ) : error("imaginary roots in quadratic")
end








###############################################################################
#
# Volume scattering functions adatpted from Volscatt version 2 by W. Verhoef
# This function passed the FT test
# This function is documented in the Utils page
#
###############################################################################
"""
    volscatt(tts::FT, tto::FT, psi::FT, ttl::FT)

Calculate interception parameters (`chi_s` and `chi_s`) and leaf reflectance multiplier (`frho`) and transmittance multiplier (`ftau`), given
- `tts` Solar zenith angle
- `tto` Viewing zenith angle
- `psi` Azimuth angle
- `ttl` Leaf inclination angle
"""
function volscatt(tts::FT, tto::FT, psi::FT, ttl::FT) where {FT}
    psi_rad = deg2rad(psi)
    cos_psi = cosd(psi)
    cos_ttl = cosd(ttl)
    sin_ttl = sind(ttl)
    cos_tts = cosd(tts)
    sin_tts = sind(tts)
    cos_tto = cosd(tto)
    sin_tto = sind(tto)
    Cs      = cos_ttl*cos_tts
    Ss      = sin_ttl*sin_tts
    Co      = cos_ttl*cos_tto
    So      = sin_ttl*sin_tto

    cosbts = FT(1.0)
    cosbto = FT(1.0)
    if (abs(Ss)>1e-6)
    	cosbts = -Cs / Ss
    end
    if (abs(So)>1e-6)
    	cosbto = -Co / So
    end

    if (abs(cosbts)<1)
    	bts = acos(cosbts)
    	ds  = Ss
    else
    	bts = FT(pi)
    	ds  = Cs
    end

    if (abs(cosbto)<1)
    	bto = acos(cosbto)
    	doo = So
    elseif(tto<90)
    	bto = FT(pi)
    	doo = Co
    else
    	bto = 0
    	doo = -Co
    end

    chi_s = 2 / FT(pi) * ( (bts - FT(pi)/2) * Cs + sin(bts) * Ss )
    chi_o = 2 / FT(pi) * ( (bto - FT(pi)/2) * Co + sin(bto) * So )

    # Computation of auxiliary azimut angles bt1, bt2, bt3 used
    # for the computation of the bidirectional scattering coefficient w

    btran1 = abs(bts - bto)
    btran2 = 2FT(pi) - bts - bto
    #btran_=pi-abs(bts+bto-FT(pi))

    if psi_rad <= btran1
    	bt1 = psi_rad
    	bt2 = btran1
    	bt3 = btran2
    else
    	bt1 = btran1
    	if psi_rad <= btran2
    		bt2 = psi_rad
    		bt3 = btran2
        else
    		bt2 = btran2
    		bt3 = psi_rad
        end
    end

    t1 = 2Cs * Co + Ss * So * cos_psi
    t2 = 0
    if bt2 > 0
    	t2 = sin(bt2) * ( 2ds * doo + Ss * So * cos(bt1) * cos(bt3) )
    end

    denom = 2 * FT(pi)^2
    frho  = ((pi-bt2) * t1 + t2) / denom
    ftau  = (-bt2     * t1 + t2) / denom

    return chi_s, abs(chi_o), frho, ftau
end








###############################################################################
#
# Functions in development
# Test and document the functions before merging them to this file
#
###############################################################################
include("MathTools_in_development.jl")




end # module
