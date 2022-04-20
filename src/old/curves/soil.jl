#=
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
=#
