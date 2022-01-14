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
            photo_set::AbstractPhotoModelParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT}
) where {FT<:AbstractFloat}
    return nothing
end




function leaf_temperature_dependence!(
            photo_set::AbstractPhotoModelParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            T::FT
) where {FT<:AbstractFloat}
    return nothing
end
