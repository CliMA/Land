###############################################################################
#
# Make g_min and g_max control
#
###############################################################################
"""
    leaf_gsw_control!(photo_set::AbstractPhotoModelParaSet{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}, ind::Int) where {FT<:AbstractFloat}
    leaf_gsw_control!(photo_set::AbstractPhotoModelParaSet{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}) where {FT<:AbstractFloat}

make sure g_sw is in its physiological range limited by diffusion, given
- `photo_set` [`C3ParaSet`] or [`C4ParaSet`] type parameter set
- `canopyi` [`CanopyLayer`](@ref) type struct
- `envir` [`AirLayer`] type struct
- `ind` Nth leaf

Note that this function is meant to use jointly with update_leaf_from_glc! when
    computing optimal stomtal conductance.
"""
function leaf_gsw_control!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            ind::Int
) where {FT<:AbstractFloat}
    # if g_sw is low than g_min
    if canopyi.g_sw[ind] < canopyi.g_min
        update_leaf_from_gsw!(photo_set, canopyi, envir, ind, canopyi.g_min);
    elseif canopyi.g_sw[ind] > canopyi.g_max
        update_leaf_from_gsw!(photo_set, canopyi, envir, ind, canopyi.g_max);
    end

    return nothing
end




function leaf_gsw_control!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT}
) where {FT<:AbstractFloat}
    # if g_sw is low than g_min
    for i in eachindex(canopyi.g_sw)
        canopyi.ps.APAR = canopyi.APAR[i];
        leaf_ETR!(photo_set, canopyi.ps);
        leaf_gsw_control!(photo_set, canopyi, envir, i);
    end

    return nothing
end
