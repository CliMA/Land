###############################################################################
#
# Temperature dependency of the photosynthetic parameters
#
###############################################################################
"""
    photo_TD_from_set(td_set::ArrheniusTD{FT}, T::FT) where {FT<:AbstractFloat}
    photo_TD_from_set(td_set::Q10TD{FT}, T::FT) where {FT<:AbstractFloat}

Make temperature correction from parameter set, given
- `td_set` [`ArrheniusTD`](@ref) type parameter set, which has a `VAL_25` field
- `T` Leaf temperature

Useful for Kc, Ko, Kpep, and ``Γ^{*}``.
"""
function photo_TD_from_set(
            td_set::ArrheniusTD{FT},
            T::FT
) where {FT<:AbstractFloat}
    return td_set.VAL_25 * temperature_correction(td_set, T)
end




function photo_TD_from_set(
            td_set::Q10TD{FT},
            T::FT
) where {FT<:AbstractFloat}
    return td_set.VAL_REF * temperature_correction(td_set, T)
end




"""
    photo_TD_from_val(
                td_set::AbstractTDParameterSet{FT},
                val::FT,
                T::FT
    ) where {FT<:AbstractFloat}

Make temperature correction from a given value, given
- `td_set` [`ArrheniusTD`](@ref) or [`ArrheniusPeakTD`](@ref) type struct
- `val` Uncorrected value at 298.15 K
- `T` Leaf temperature

Useful for Vcmax, Vomax, Vpmax, Jmax, and Respiration.
"""
function photo_TD_from_val(
            td_set::AbstractTDParameterSet{FT},
            val::FT,
            T::FT
) where {FT<:AbstractFloat}
    return val * temperature_correction(td_set, T)
end








###############################################################################
#
# Functions to update the TD individually
#
###############################################################################
"""
    leaf_jmax!(td_set::AbstractTDParameterSet{FT},
               leaf::Leaf{FT}
    ) where {FT<:AbstractFloat}

Update maximal electron transport rate at leaf temperature, given
- `td_set` [`AbstractTDParameterSet`](@ref) type TD parameter set
- `leaf` [`Leaf`](@ref) type struct
"""
function leaf_jmax!(
            td_set::AbstractTDParameterSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    leaf.Jmax = photo_TD_from_val(td_set, leaf.Jmax25, leaf.T);

    return nothing
end




"""
    leaf_kc!(td_set::ArrheniusTD{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}

Update Kc at leaf temperature, given
- `td_set` [`ArrheniusTD`](@ref) type TD parameter set
- `leaf` [`Leaf`](@ref) type struct
"""
function leaf_kc!(
            td_set::ArrheniusTD{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    leaf.Kc = photo_TD_from_set(td_set, leaf.T);

    return nothing
end




"""
    leaf_km!(photo_set::C3ParaSet{FT},
             leaf::Leaf{FT},
             envir::AirLayer{FT}
    ) where {FT<:AbstractFloat}

Update Ko at leaf temperature, given
- `photo_set` [`C3ParaSet`](@ref) type photosynthesis parameter set
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
"""
function leaf_km!(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT}
) where {FT<:AbstractFloat}
    leaf.Km = leaf.Kc * (1 + envir.p_O₂/leaf.Ko);

    return nothing
end




"""
    leaf_ko!(td_set::ArrheniusTD{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}

Update Ko at leaf temperature, given
- `td_set` [`ArrheniusTD`](@ref) type TD parameter set
- `leaf` [`Leaf`](@ref) type struct
"""
function leaf_ko!(
            td_set::ArrheniusTD{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    leaf.Ko = photo_TD_from_set(td_set, leaf.T);

    return nothing
end




"""
    leaf_kpep!(td_set::ArrheniusTD{FT},
               leaf::Leaf{FT}
    ) where {FT<:AbstractFloat}

Update Kpep at leaf temperature, given
- `td_set` [`ArrheniusTD`](@ref) type TD parameter set
- `leaf` [`Leaf`](@ref) type struct
"""
function leaf_kpep!(
            td_set::ArrheniusTD{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    leaf.Kpep = photo_TD_from_set(td_set, leaf.T);

    return nothing
end




"""
    leaf_rd!(td_set::AbstractTDParameterSet{FT},
             leaf::Leaf{FT}
    ) where {FT<:AbstractFloat}

Update leaf dark respiration rate at leaf temperature, given
- `td_set` [`AbstractTDParameterSet`](@ref) type TD parameter set
- `leaf` [`Leaf`](@ref) type struct
"""
function leaf_rd!(
            td_set::AbstractTDParameterSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    leaf.Rd = photo_TD_from_val(td_set, leaf.Rd25, leaf.T);

    return nothing
end




"""
    leaf_vcmax!(td_set::AbstractTDParameterSet{FT},
                leaf::Leaf{FT}
    ) where {FT<:AbstractFloat}

Update leaf maximal carboxylation rate at leaf temperature, given
- `td_set` [`AbstractTDParameterSet`](@ref) type TD parameter set
- `leaf` [`Leaf`](@ref) type struct
"""
function leaf_vcmax!(
            td_set::AbstractTDParameterSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    leaf.Vcmax = photo_TD_from_val(td_set, leaf.Vcmax25, leaf.T);

    return nothing
end




"""
    leaf_vpmax!(td_set::AbstractTDParameterSet{FT},
                leaf::Leaf{FT}
    ) where {FT<:AbstractFloat}

Update leaf maximal PEP carboxylation rate at leaf temperature, given
- `td_set` [`AbstractTDParameterSet`](@ref) type TD parameter set
- `leaf` [`Leaf`](@ref) type struct
"""
function leaf_vpmax!(
            td_set::AbstractTDParameterSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    leaf.Vpmax = photo_TD_from_val(td_set, leaf.Vpmax25, leaf.T);

    return nothing
end




"""
    leaf_Γstar!(td_set::ArrheniusTD{FT},
                leaf::Leaf{FT}
    ) where {FT<:AbstractFloat}

Update ``Γ^{*}`` at leaf temperature, given
- `td_set` [`ArrheniusTD`](@ref) type TD parameter set
- `leaf` [`Leaf`](@ref) type struct
"""
function leaf_Γstar!(
            td_set::ArrheniusTD{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    leaf.Γ_star = photo_TD_from_set(td_set, leaf.T);

    return nothing
end








###############################################################################
#
# Calculate photosynthesis using Leaf
#
###############################################################################
"""
    leaf_temperature_dependence!(
                photo_set::AbstractPhotoModelParaSet{FT},
                leaf::Leaf{FT},
                envir::AirLayer{FT}
    ) where {FT<:AbstractFloat}
    leaf_temperature_dependence!(
                photo_set::AbstractPhotoModelParaSet{FT},
                leaf::Leaf{FT},
                envir::AirLayer{FT},
                T::FT
    ) where {FT<:AbstractFloat}

Update the temperature dependent photosynthesis only, given
- `photo_set` [`AbstractPhotoModelParaSet`](@ref) type parameter set
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
- `T` Given leaf temperature
"""
function leaf_temperature_dependence!(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT}
) where {FT<:AbstractFloat}
    if leaf.T_old != leaf.T
        leaf.T_old = leaf.T;
        leaf.p_sat = saturation_vapor_pressure(leaf.T);
        leaf_rd!(photo_set.ReT, leaf);
        leaf_vcmax!(photo_set.VcT, leaf);
        leaf_jmax!(photo_set.JT , leaf);
        leaf_kc!(photo_set.KcT, leaf);
        leaf_ko!(photo_set.KoT, leaf);
        leaf_km!(photo_set, leaf, envir);
        leaf_Γstar!(photo_set.ΓsT, leaf);
    end

    return nothing
end




function leaf_temperature_dependence!(
            photo_set::C4ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT}
) where {FT<:AbstractFloat}
    if leaf.T_old != leaf.T
        leaf.T_old = leaf.T
        leaf.p_sat = saturation_vapor_pressure(leaf.T);
        leaf_rd!(photo_set.ReT, leaf);
        leaf_vcmax!(photo_set.VcT, leaf);
        leaf_vpmax!(photo_set.VpT, leaf);
        leaf_kpep!(photo_set.KpT, leaf);
    end

    return nothing
end




function leaf_temperature_dependence!(
            photo_set::AbstractPhotoModelParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            T::FT
) where {FT<:AbstractFloat}
    leaf.T = T;
    leaf_temperature_dependence!(photo_set, leaf, envir);

    return nothing
end
