###############################################################################
#
# Update CanopyLayer struct from stomatal conductance
#
###############################################################################
"""
    update_leaf_from_glc!(photo_set::AbstractPhotoModelParaSet{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}, ind::Int, glc::FT) where {FT<:AbstractFloat}
    update_leaf_from_glc!(photo_set::AbstractPhotoModelParaSet{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}) where {FT<:AbstractFloat}

Update Nth leaf photosynthesis, given
- `photo_set` [`C3ParaSet`] or [`C4ParaSet`] type parameter set
- `canopyi` [`CanopyLayer`](@ref) type struct
- `envir` [`AirLayer`] type struct
- `ind` Nth leaf
- `glc` Given leaf diffusive conductance

Note that this function does not make the gsw control, so it is not guaranteed
    that gsw is within the physiological range. Thus, gsw control should be
    made outside this function. This function is supposed to be used in the
    optimal stomatl conductance models only, because optimal conductance can be
    outside the physiological stomatal conductance range. Thus, using this
    function for other purposes need to be cautious. In this case, it is
    recommended to use `update_leaf_from_gsw!`.
"""
function update_leaf_from_glc!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            ind::Int,
            glc::FT
) where {FT<:AbstractFloat}
    # update the conductances
    canopyi.g_lc[ind] = glc;
    update_leaf_from_glc!(photo_set, canopyi, envir, ind);

    return nothing
end




function update_leaf_from_glc!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            ind::Int
) where {FT<:AbstractFloat}
    # update the conductances
    canopyi.g_sc[ind] = 1 / ( 1 / canopyi.g_lc[ind] -
                              1 / canopyi.g_m[ind] -
                              1 / canopyi.g_bc[ind] );
    canopyi.g_sw[ind] = canopyi.g_sc[ind] * FT(1.6);
    canopyi.g_lw[ind] = 1 / ( 1 / canopyi.g_sw[ind] +
                              1 / canopyi.g_bw[ind] );

    # update the photosynthetic rates
    if canopyi.g_lc[ind] != canopyi.ps.g_lc
        leaf_photo_from_glc!(photo_set, canopyi.ps, envir, canopyi.g_lc[ind]);
        #
        #
        #
        #
        # be careful that this one might have memory allocation
        # make some changes on Photosynthesis.jl to avoid memory allocation
        # such as leaf_fluorescence!(photo_set, canopyi.ps);
        #
        #
        #
        #
        leaf_fluorescence!(photo_set.Flu, canopyi.ps);
    end
    canopyi.Ac[ind] = canopyi.ps.Ac;
    canopyi.Aj[ind] = canopyi.ps.Aj;
    canopyi.Ap[ind] = canopyi.ps.Ap;
    canopyi.Ag[ind] = canopyi.ps.Ag;
    canopyi.An[ind] = canopyi.ps.An;
    canopyi.ϕs[ind] = canopyi.ps.ϕs;

    # update the pressures
    canopyi.p_i[ind] = canopyi.ps.p_i;
    canopyi.p_s[ind] = canopyi.ps.p_s;

    return nothing
end




function update_leaf_from_glc!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT}
) where {FT<:AbstractFloat}
    # update the conductances for each "leaf"
    for i in eachindex(canopyi.g_lc)
        canopyi.ps.APAR = canopyi.APAR[i];
        leaf_ETR!(photo_set, canopyi.ps);
        update_leaf_from_glc!(photo_set, canopyi, envir, i, canopyi.g_lc[i]);
    end

    return nothing
end




"""
    update_leaf_from_gsw!(photo_set::AbstractPhotoModelParaSet{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}, ind::Int, gsw::FT) where {FT<:AbstractFloat}
    update_leaf_from_gsw!(photo_set::AbstractPhotoModelParaSet{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}, ind::Int) where {FT<:AbstractFloat}
    update_leaf_from_gsw!(photo_set::AbstractPhotoModelParaSet{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}) where {FT<:AbstractFloat}

Update Nth leaf photosynthesis, given
- `photo_set` [`C3ParaSet`] or [`C4ParaSet`] type parameter set
- `canopyi` [`CanopyLayer`](@ref) type struct
- `envir` [`AirLayer`] type struct
- `ind` Nth leaf
- `gsw` Given stomatal conductance to H₂O

Note that this function makes the gsw control so that gsw is within the
    physiological range.
"""
function update_leaf_from_gsw!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            ind::Int,
            gsw::FT
) where {FT<:AbstractFloat}
    # update the conductances
    canopyi.g_sw[ind] = gsw;
    update_leaf_from_gsw!(photo_set, canopyi, envir, ind);

    return nothing
end




function update_leaf_from_gsw!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            ind::Int
) where {FT<:AbstractFloat}
    @unpack g_max, g_min = canopyi;

    # update the conductances
    if canopyi.g_sw[ind] > g_max
        canopyi.g_sw[ind] = g_max;
    end

    if canopyi.g_sw[ind] < g_min
        canopyi.g_sw[ind] = g_min;
    end

    canopyi.g_lw[ind] = 1 / ( 1 / canopyi.g_sw[ind] +
                              1 / canopyi.g_bw[ind] );
    canopyi.g_sc[ind] = canopyi.g_sw[ind] / FT(1.6);
    canopyi.g_lc[ind] = 1 / ( FT(1.6) / canopyi.g_sw[ind] +
                              1 / canopyi.g_m[ind] +
                              1 / canopyi.g_bc[ind] );

    # update the photosynthetic rates
    leaf_photo_from_glc!(photo_set, canopyi.ps, envir, canopyi.g_lc[ind]);
    canopyi.Ac[ind] = canopyi.ps.Ac;
    canopyi.Aj[ind] = canopyi.ps.Aj;
    canopyi.Ap[ind] = canopyi.ps.Ap;
    canopyi.Ag[ind] = canopyi.ps.Ag;
    canopyi.An[ind] = canopyi.ps.An;

    # update the pressures
    canopyi.p_i[ind] = canopyi.ps.p_i;
    canopyi.p_s[ind] = canopyi.ps.p_s;

    return nothing
end




function update_leaf_from_gsw!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT}
            ) where {FT<:AbstractFloat}
    # update the conductances for each "leaf"
    for i in eachindex(canopyi.g_sw)
        canopyi.ps.APAR = canopyi.APAR[i];
        leaf_ETR!(photo_set, canopyi.ps);
        update_leaf_from_gsw!(photo_set, canopyi, envir, i);
    end

    return nothing
end
