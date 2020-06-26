module MathTools

using IncGammaBeta
using Parameters

using ..LandParameters

const AVOGADRO    = LandParameters.AVOGADRO
const H_PLANCK    = LandParameters.H_PLANCK
const LIGHT_SPEED = LandParameters.LIGHT_SPEED

export dladgen,
       e2phot,
       fast∫,
       lower_quadratic,
       psofunction,
       quadratic,
       volscatt,
       weibull_k_ratio




###############################################################################
#
# Convert energy to photons
#
###############################################################################
const _FAC = 1e-9 / (H_PLANCK * LIGHT_SPEED * AVOGADRO)

"""
    e2phot(λ::Array{FT},E::Array{FT})

Calculates the number of moles of photons, given
- `λ` An array of wave length in `[nm]`, converted to `[m]` by _FAC
- `E` Joules of energy
"""
function e2phot(
            λ::Array{FT},
            E::Array{FT}
            ) where {FT<:AbstractFloat}
    return (E .* λ) .* FT(_FAC)
end








###############################################################################
#
# Integral functions
#
###############################################################################
"""
    fast∫(dx::Vector{FT}, f::Vector{FT})

A fast way of integrating functions, given
- `dx` Delta x for each x
- `f` f(x) for each x
"""
function fast∫(
            dx::Vector{FT},
            f::Vector{FT}
            ) where {FT<:AbstractFloat}
    if length(dx) == length(f)
        return sum( f .* dx )
    else
        N = length(f);
        result::FT = 0.0;
        @inbounds for i=1:N
            result += f[i] * dx[i];
        end
        return result
    end
end








###############################################################################
#
# Quadratic function solver
#
###############################################################################
"""
    lower_quadratic(a::FT, b::FT, c::FT)

Return the lower quadratic solution or NaN, given
- `a` Parameter in `a*x^2 + b*x + c = 0`
- `b` Parameter in `a*x^2 + b*x + c = 0`
- `c` Parameter in `a*x^2 + b*x + c = 0`
"""
function lower_quadratic(
            a::FT,
            b::FT,
            c::FT
            ) where {FT<:AbstractFloat}
    discr = b^2 - 4*a*c;
    discr >= 0 ? (-b - sqrt(discr))/2a : FT(NaN)
end




"""
    quadratic(a::FT, b::FT, c::FT)

Return a list of quadratic function solutions or an error message, given
- `a` Parameter in `a*x^2 + b*x + c = 0`
- `b` Parameter in `a*x^2 + b*x + c = 0`
- `c` Parameter in `a*x^2 + b*x + c = 0`
"""
function quadratic(
            a::FT,
            b::FT,
            c::FT
            ) where {FT<:AbstractFloat}
    discr = b^2 - 4*a*c;
    discr >= 0 ? ( (-b + sqrt(discr))/2a, (-b - sqrt(discr))/2a ) :
                 error("imaginary roots in quadratic")
end








###############################################################################
#
# Volume scattering functions adatpted from Volscatt version 2 by W. Verhoef
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
function volscatt(
            tts::FT,
            tto::FT,
            psi::FT,
            ttl::FT
            ) where {FT<:AbstractFloat}
    psi_rad = deg2rad(psi);
    cos_psi = cosd(psi);
    cos_ttl = cosd(ttl);
    sin_ttl = sind(ttl);
    cos_tts = cosd(tts);
    sin_tts = sind(tts);
    cos_tto = cosd(tto);
    sin_tto = sind(tto);
    Cs      = cos_ttl*cos_tts;
    Ss      = sin_ttl*sin_tts;
    Co      = cos_ttl*cos_tto;
    So      = sin_ttl*sin_tto;

    cosbts = FT(1.0);
    cosbto = FT(1.0);
    if (abs(Ss)>1e-6)
    	cosbts = -Cs / Ss;
    end
    if (abs(So)>1e-6)
    	cosbto = -Co / So;
    end

    if (abs(cosbts)<1)
    	bts = acos(cosbts);
    	ds  = Ss;
    else
    	bts = FT(pi);
    	ds  = Cs;
    end

    if (abs(cosbto)<1)
    	bto = acos(cosbto);
    	doo = So;
    elseif(tto<90)
    	bto = FT(pi);
    	doo = Co;
    else
    	bto = 0;
    	doo = -Co;
    end

    chi_s = 2 / FT(pi) * ( (bts - FT(pi)/2) * Cs + sin(bts) * Ss );
    chi_o = 2 / FT(pi) * ( (bto - FT(pi)/2) * Co + sin(bto) * So );

    # Computation of auxiliary azimut angles bt1, bt2, bt3 used
    # for the computation of the bidirectional scattering coefficient w

    btran1 = abs(bts - bto);
    btran2 = 2FT(pi) - bts - bto;
    #btran_=pi-abs(bts+bto-FT(pi))

    if psi_rad <= btran1
    	bt1 = psi_rad;
    	bt2 = btran1;
    	bt3 = btran2;
    else
    	bt1 = btran1;
    	if psi_rad <= btran2
    		bt2 = psi_rad;
    		bt3 = btran2;
        else
    		bt2 = btran2;
    		bt3 = psi_rad;
        end
    end

    t1 = 2Cs * Co + Ss * So * cos_psi;
    t2 = 0;
    if bt2 > 0
    	t2 = sin(bt2) * ( 2ds * doo + Ss * So * cos(bt1) * cos(bt3) );
    end

    denom = 2 * FT(pi)^2;
    frho  = ((pi-bt2) * t1 + t2) / denom;
    ftau  = (-bt2     * t1 + t2) / denom;

    return chi_s, abs(chi_o), frho, ftau
end








###############################################################################
#
# Weibull functions
#
###############################################################################
"""
    weibull_k_ratio(b::FT, c::FT, p_25::FT, vis::FT)

Returns the relative hydraulic conductance, given
- `b` Weibull B in `[MPa]`
- `c` Weibull C
- `p_25` Equivalent xylem pressure at 298.15 K in `[MPa]`
- `vis` Relative viscosity
"""
function weibull_k_ratio(
            b::FT,
            c::FT,
            p_25::FT,
            vis::FT
            ) where {FT<:AbstractFloat}
    kr = max( FT(1e-4), exp( -1 * (-p_25/b) ^ c ) / vis );

    return kr
end








###############################################################################
#
# Functions in development
# Test and document the functions before merging them to this file
#
###############################################################################
include("MathTools_in_development.jl")




end # module
