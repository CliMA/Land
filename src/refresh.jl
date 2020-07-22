###############################################################################
#
# Update CanopyLayer struct when temperature changes
#
###############################################################################
"""
    update_leaf_TP!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            envir::AirLayer{FT}) where {FT<:AbstractFloat}

Update leaf physiological parameters if temperature or pressure changes in the
daytime, given
- `photo_set` [`C3ParaSet`] or [`C4ParaSet`] type parameter set
- `canopyi` [`CanopyLayer`](@ref) type struct
- `hs` Leaf hydraulic system
- `envir` [`AirLayer`] type struct
"""
function update_leaf_TP!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            envir::AirLayer{FT}
            ) where {FT<:AbstractFloat}
    # unpack required variables
    @unpack g_max25, g_min25, p_ups, p_old, T, T_old = canopyi;

    # if T changes, update TD and ec, then T_old
    if T != T_old
        canopyi.g_max    = g_max25 * relative_diffusive_coefficient(T);
        canopyi.g_min    = g_min25 * relative_diffusive_coefficient(T);
        canopyi.LV       = latent_heat_vapor(T) * 1000 / FT(MOLMASS_WATER);
        canopyi.ps.T     = canopyi.T;
        leaf_temperature_dependence!(photo_set, canopyi.ps, envir);
        canopyi.p_sat    = canopyi.ps.p_sat;
        hs.f_st         = relative_surface_tension(T);
        hs.f_vis        = relative_viscosity(T);
        hs.p_ups        = canopyi.p_ups;
        canopyi.ec       = leaf_e_crit(hs, canopyi.ec);
        canopyi.T_old    = canopyi.T;
        canopyi.p_old    = canopyi.p_ups;
    # if only p_ups changes, update ec and p_old
    elseif p_ups != p_old
        hs.p_ups = canopyi.p_ups;
        canopyi.ec       = leaf_e_crit(hs, canopyi.ec);
        canopyi.p_old    = canopyi.p_ups;
    end

    return nothing
end




"""
    update_leaf_AK!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            envir::AirLayer{FT}) where {FT<:AbstractFloat}

Update leaf maximal A and K for Sperry model, given
- `photo_set` [`C3ParaSet`] or [`C4ParaSet`] type parameter set
- `canopyi` [`CanopyLayer`](@ref) type struct
- `envir` [`AirLayer`] type struct
"""
function update_leaf_AK!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            envir::AirLayer{FT}
            ) where {FT<:AbstractFloat}
    # unpack required variables
    @unpack APAR, ec, g_bc, g_bw, g_m, g_max, g_min, p_sat = canopyi;
    @unpack p_atm, p_H₂O = envir;

    # calculate the physiological maximal g_sw
    _g_crit = ec / (p_sat - p_H₂O) * p_atm;
    _g_sw   = 1 ./ max.(1/_g_crit .- 1 ./ g_bw, FT(1e-3));
    _g_sw  .= min.(_g_sw, g_max);
    _g_sw  .= max.(_g_sw, g_min);
    _g_lcs  = 1 ./ (1 ./ g_bc .+ FT(1.6) ./ _g_sw .+ 1 ./ g_m);

    # update the limited rates
    _Jps = leaf_ETR_pot_APAR(canopyi.ps, APAR);
    _Js  = leaf_ETR_Jps(photo_set, canopyi.ps, _Jps);
    _Anc = rubisco_limited_an_glc(photo_set, canopyi.ps, envir, _g_lcs);
    _Anj = light_limited_an_glc(photo_set, canopyi.ps, envir, _g_lcs, _Js);
    _Anp = product_limited_an_glc(photo_set, canopyi.ps, envir, _g_lcs);
    _Ams = min.(_Anc, _Anj, _Anp);

    # update a_max and kr_max
    canopyi.a_max  .= _Ams;
    canopyi.kr_max  = hs.k_history[end];

    return nothing
end








###############################################################################
#
# Update CanopyLayer struct from stomatal conductance
#
###############################################################################
"""
    update_leaf_from_glc!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            ind::Int,
            glc::FT)

Update Nth leaf photosynthesis, given
- `photo_set` [`C3ParaSet`] or [`C4ParaSet`] type parameter set
- `canopyi` [`CanopyLayer`](@ref) type struct
- `envir` [`AirLayer`] type struct
- `ind` Nth leaf
- `glc` Given leaf diffusive conductance
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
    canopyi.g_sc[ind] = 1 / ( 1 / glc -
                             1 / canopyi.g_m[ind] -
                             1 / canopyi.g_bc[ind] );
    canopyi.g_sw[ind] = canopyi.g_sc[ind] * FT(1.6);
    canopyi.g_lw[ind] = 1 / ( 1 / canopyi.g_sw[ind] +
                             1 / canopyi.g_bw[ind] );

    # update the photosynthetic rates
    if glc != canopyi.ps.g_lc
        leaf_photo_from_glc!(photo_set, canopyi.ps, envir, glc);
    end
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




function update_leaf_from_glc!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT}
            ) where {FT<:AbstractFloat}
    # update the conductances for each "leaf"
    for i in eachindex(canopyi.g_lc)
        update_leaf_from_glc!(photo_set, canopyi, envir, i, canopyi.g_lc[i]);
    end

    return nothing
end




"""
    update_leaf_from_gsw!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            ind::Int,
            glc::FT)

Update Nth leaf photosynthesis, given
- `photo_set` [`C3ParaSet`] or [`C4ParaSet`] type parameter set
- `canopyi` [`CanopyLayer`](@ref) type struct
- `envir` [`AirLayer`] type struct
- `ind` Nth leaf
- `gsw` Given stomatal conductance to H₂O
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
    canopyi.g_lw[ind] = 1 / ( 1 / gsw +
                             1 / canopyi.g_bw[ind] );
    canopyi.g_sc[ind] = gsw / FT(1.6);
    canopyi.g_lc[ind] = 1 / ( FT(1.6) / gsw +
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
        update_leaf_from_gsw!(photo_set, canopyi, envir, i, canopyi.g_sw[i]);
    end

    return nothing
end








###############################################################################
#
# Make g_min and g_max control
#
###############################################################################
"""
    leaf_gsw_control!(
            photo_set::AbstractPhotoModelParaSet,
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            ind::Int)

make sure g_sw is in its physiological range limited by diffusion, given
- `photo_set` [`C3ParaSet`] or [`C4ParaSet`] type parameter set
- `canopyi` [`CanopyLayer`](@ref) type struct
- `envir` [`AirLayer`] type struct
"""
function leaf_gsw_control!(
            photo_set::AbstractPhotoModelParaSet,
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
            photo_set::AbstractPhotoModelParaSet,
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
