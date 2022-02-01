###############################################################################
#
# Equilibrium pressure from volume
#
###############################################################################
"""
    p_from_volume(
                pv::AbstractCapacity{FT},
                rvol::FT,
                T::FT
    ) where {FT<:AbstractFloat}

Calculate equilibrium pressure from relative volume, given
- `pv` [`AbstractCapacity`](@ref) type of struct
"""
function p_from_volume(
            pv::PVCurveLinear{FT},
            rvol::FT,
            T::FT
) where {FT<:AbstractFloat}
    return (rvol - 1) / pv.slope
end




function p_from_volume(
            pv::PVCurveSegmented{FT},
            rvol::FT,
            T::FT
) where {FT<:AbstractFloat}
    @unpack c_all, RWC_apo, RWC_TLP, ϵ_bulk = pv;

    if rvol > RWC_TLP
        return -c_all * GAS_R(FT) * T / (rvol - RWC_apo) * FT(1e-6) +
                ϵ_bulk * (rvol - RWC_TLP)
    elseif rvol > RWC_apo
        return -c_all * GAS_R(FT) * T / (rvol - RWC_apo) * FT(1e-6)
    else
        return FT(-100)
    end
end
