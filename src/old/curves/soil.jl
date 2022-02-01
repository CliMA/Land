###############################################################################
#
# Calculate soil water content from pressure
#
###############################################################################
"""
    soil_erwc(sh::AbstractSoilVC{FT}, p_25::FT) where {FT<:AbstractFloat}

Returns the Effective relative water content of the soil
    ``\\frac{Θ - Θr}{Θs - Θr}``, given
- `sh` [`BrooksCorey`](@ref) or [`VanGenuchten`](@ref) type soil hydraulics
- `p_25` Matrix water potential equivalent to 25 degree C, with surface tension
correction
"""
function soil_erwc(
            sh::BrooksCorey{FT},
            p_25::FT
) where {FT<:AbstractFloat}
    if p_25 < 0
        @unpack b, ϕs = sh;

        return (-ϕs/p_25) ^ (1/b)
    else
        return FT(1)
    end
end




function soil_erwc(
            sh::VanGenuchten{FT},
            p_25::FT
) where {FT<:AbstractFloat}
    if p_25 < 0
        @unpack m, n, α = sh;

        return ( 1 / ( 1 + (-p_25 * α) ^ n ) ) ^ m
    else
        return FT(1)
    end
end




"""
    soil_rwc(sh::AbstractSoilVC{FT}, p_25::FT) where {FT<:AbstractFloat}

Returns the relative soil water content, given
- `sh` [`BrooksCorey`](@ref) or [`VanGenuchten`](@ref) type soil hydraulics
- `p_25` Matrix water potential equivalent to 25 degree C, with surface tension
correction
"""
function soil_rwc(
            sh::BrooksCorey{FT},
            p_25::FT
) where {FT<:AbstractFloat}
    @unpack Θs = sh;

    return soil_swc(sh, p_25) / Θs
end




function soil_rwc(
            sh::VanGenuchten{FT},
            p_25::FT
) where {FT<:AbstractFloat}
    @unpack Θs = sh;

    return soil_swc(sh, p_25) / Θs
end




"""
    soil_swc(sh::AbstractSoilVC{FT}, p_25::FT) where {FT<:AbstractFloat}

Returns the soil water content, given
- `sh` [`BrooksCorey`](@ref) or [`VanGenuchten`](@ref) type soil hydraulics
- `p_25` Matrix water potential equivalent to 25 degree C, with surface tension
correction
"""
function soil_swc(
            sh::BrooksCorey{FT},
            p_25::FT
) where {FT<:AbstractFloat}
    @unpack Θr, Θs = sh;

    return soil_erwc(sh, p_25) * (Θs - Θr) + Θr
end




function soil_swc(
            sh::VanGenuchten{FT},
            p_25::FT
) where {FT<:AbstractFloat}
    @unpack Θr, Θs = sh;

    return soil_erwc(sh, p_25) * (Θs - Θr) + Θr
end








###############################################################################
#
# Calculate soil k from water content
#
###############################################################################
"""
    soil_k_ratio_erwc(
                sh::AbstractSoilVC{FT},
                erwc::FT
    ) where {FT<:AbstractFloat}

Return the soil k relative to maximal k, given
- `sh` [`BrooksCorey`](@ref) or [`VanGenuchten`](@ref) type soil hydraulics
- `erwc` Effective relative soil water content (``\\frac{Θs - Θr}{Θs - Θ}``)
"""
function soil_k_ratio_erwc(
            sh::BrooksCorey{FT},
            erwc::FT
) where {FT<:AbstractFloat}
    k_ratio = erwc ^ (2*sh.b+3);

    return max(k_ratio, FT(1e-20))
end




function soil_k_ratio_erwc(
            sh::VanGenuchten{FT},
            erwc::FT
) where {FT<:AbstractFloat}
    @unpack m = sh;
    k_ratio = sqrt(erwc) * (1 - (1 - erwc^(1/m)) ^ m)^2;

    return max(k_ratio, FT(1e-20))
end




"""
    soil_k_ratio_rwc(sh::AbstractSoilVC{FT}, rwc::FT) where {FT<:AbstractFloat}

Return the soil k relative to maximal k, given
- `sh` [`BrooksCorey`](@ref) or [`VanGenuchten`](@ref) type soil hydraulics
- `rwc` Relative soil water content
"""
function soil_k_ratio_rwc(
            sh::BrooksCorey{FT},
            rwc::FT
) where {FT<:AbstractFloat}
    @unpack Θr, Θs = sh;

    erwc = min(1, max(0, (rwc * Θs - Θr) / (Θs - Θr)));

    return soil_k_ratio_erwc(sh, erwc)
end




function soil_k_ratio_rwc(
            sh::VanGenuchten{FT},
            rwc::FT
) where {FT<:AbstractFloat}
    @unpack Θr, Θs = sh;

    erwc = min(1, max(0, (rwc * Θs - Θr) / (Θs - Θr)));

    return soil_k_ratio_erwc(sh, erwc)
end




"""
    soil_k_ratio_swc(sh::AbstractSoilVC{FT}, swc::FT) where {FT<:AbstractFloat}

Return the soil k relative to maximal k, given
- `sh` [`BrooksCorey`](@ref) or [`VanGenuchten`](@ref) type soil hydraulics
- `swc` Relative soil water content
"""
function soil_k_ratio_swc(
            sh::BrooksCorey{FT},
            swc::FT
) where {FT<:AbstractFloat}
    @unpack Θr, Θs = sh;

    erwc = min(1, max(0, (swc - Θr) / (Θs - Θr)));

    return soil_k_ratio_erwc(sh, erwc)
end




function soil_k_ratio_swc(
            sh::VanGenuchten{FT},
            swc::FT
) where {FT<:AbstractFloat}
    @unpack Θr, Θs = sh;

    erwc = min(1, max(0, (swc - Θr) / (Θs - Θr)));

    return soil_k_ratio_erwc(sh, erwc)
end








###############################################################################
#
# Calculate soil k from pressure
#
###############################################################################
"""
    soil_k_ratio_p25(
                sh::AbstractSoilVC{FT},
                p_25::FT
    ) where {FT<:AbstractFloat}

Return the soil k relative to maximal k, given
- `sh` [`AbstractSoilVC`](@ref) type soil hydraulics
- `p_25` Matrix water potential equivalent to 25 degree C, with surface tension
"""
function soil_k_ratio_p25(
            sh::AbstractSoilVC{FT},
            p_25::FT
) where {FT<:AbstractFloat}
    erwc = soil_erwc(sh, p_25);

    return soil_k_ratio_erwc(sh, erwc)
end








###############################################################################
#
# Calculate soil pressure from water content
#
###############################################################################
"""
    soil_p_25_erwc(sh::AbstractSoilVC{FT}, erwc::FT) where {FT<:AbstractFloat}

Returns the Relative water content of the soil, given
- `sh` [`BrooksCorey`](@ref) or [`VanGenuchten`](@ref) type soil hydraulics
- `erwc` Effectibe relative soil water content (``\\frac{Θs - Θr}{Θs - Θ}``)

Note that this function returns the matrix potential of water, not water
potential. Also, the potential is that at 25 degree C, not yet been corrected
for temperature effect on surface tension.
"""
function soil_p_25_erwc(
            sh::BrooksCorey{FT},
            erwc::FT
) where {FT<:AbstractFloat}
    if erwc < 1
        @unpack b, ϕs = sh;

        return -ϕs / (erwc ^ b)
    else
        return FT(0)
    end
end




function soil_p_25_erwc(
            sh::VanGenuchten{FT},
            erwc::FT
) where {FT<:AbstractFloat}
    if erwc < 1
        @unpack m, n, α = sh;

        return -1 * (erwc ^ (-1/m) - 1) ^ (1/n) / α
    else
        return FT(0)
    end
end




"""
    soil_p_25_rwc(sh::AbstractSoilVC{FT}, rwc::FT) where {FT<:AbstractFloat}

Returns the Relative water content of the soil, given
- `sh` [`BrooksCorey`](@ref) or [`VanGenuchten`](@ref) type soil hydraulics
- `rwc` Relative soil water content

Note that this function returns the matrix potential of water, not water
potential. Also, the potential is that at 25 degree C, not yet been corrected
for temperature effect on surface tension.
"""
function soil_p_25_rwc(
            sh::BrooksCorey{FT},
            rwc::FT
) where {FT<:AbstractFloat}
    @unpack Θr, Θs = sh;

    erwc = min(1, max(0, (rwc * Θs - Θr) / (Θs - Θr)));

    return soil_p_25_erwc(sh, erwc)
end




function soil_p_25_rwc(
            sh::VanGenuchten{FT},
            rwc::FT
) where {FT<:AbstractFloat}
    @unpack Θr, Θs = sh;

    erwc = min(1, max(0, (rwc * Θs - Θr) / (Θs - Θr)));

    return soil_p_25_erwc(sh, erwc)
end




"""
    soil_p_25_swc(sh::AbstractSoilVC{FT}, rwc::FT) where {FT<:AbstractFloat}

Returns the Relative water content of the soil, given
- `sh` [`BrooksCorey`](@ref) or [`VanGenuchten`](@ref) type soil hydraulics
- `swc` Soil water content

Note that this function returns the matrix potential of water, not water
potential. Also, the potential is that at 25 degree C, not yet been corrected
for temperature effect on surface tension.
"""
function soil_p_25_swc(
            sh::BrooksCorey{FT},
            swc::FT
) where {FT<:AbstractFloat}
    @unpack Θr, Θs = sh;

    erwc = min(1, max(0, (swc - Θr) / (Θs - Θr)));

    return soil_p_25_erwc(sh, erwc)
end




function soil_p_25_swc(
            sh::VanGenuchten{FT},
            swc::FT
) where {FT<:AbstractFloat}
    @unpack Θr, Θs = sh;

    erwc = min(1, max(0, (swc - Θr) / (Θs - Θr)));

    return soil_p_25_erwc(sh, erwc)
end








###############################################################################
#
# Fit Soil VC
#
###############################################################################
"""
    fit_soil_VC!(
                vc_vG::VanGenuchten{FT},
                vc_BC::BrooksCorey{FT}
    ) where {FT<:AbstractFloat}

Update BrooksCorey setup from known VanGenuchten parameters, given
- `vc_vG` [`VanGenuchten`](@ref) type soil VC
- `vc_BC` [`BrooksCorey`](@ref) type soil VC
"""
function fit_soil_VC!(
            vc_vG::VanGenuchten{FT},
            vc_BC::BrooksCorey{FT}
) where {FT<:AbstractFloat}
    # generate data to fit
    _Θ    = range(vc_vG.Θr+FT(1e-2); stop=vc_vG.Θs-FT(1e-2), length=30);
    _P_vG = -1 .* soil_p_25_swc.([vc_vG], _Θ);
    _P_BC = -1 .* soil_p_25_swc.([vc_BC], _Θ);

    f(x) = (
        vc_BC.b  = x[1];
        vc_BC.ϕs = x[2];
        _P_BC   .= -1 .* soil_p_25_swc.([vc_BC], _Θ);
        _diff    = sum( (log.(_P_BC) .- log.(_P_vG)) .^ 2 );
        return -_diff
    );

    st = SolutionToleranceND{FT}([1e-3, 1e-6], 30);
    ms = ReduceStepMethodND{FT}(x_mins = FT[1e-3, 1e-6],
                                x_maxs = FT[ 100, 1000],
                                x_inis = [vc_BC.b, vc_BC.ϕs],
                                Δ_inis = FT[0.1, 1e-3]);
    bc = find_peak(f, ms, st);
    vc_BC.b  = bc[1];
    vc_BC.ϕs = bc[2];

    return nothing
end
