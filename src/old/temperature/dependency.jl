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
        ##  leaf_rd!(photo_set.ReT, leaf);
        ##  leaf_vcmax!(photo_set.VcT, leaf);
        ##  leaf_vpmax!(photo_set.VpT, leaf);
        ##  leaf_kpep!(photo_set.KpT, leaf);
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
