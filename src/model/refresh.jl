###############################################################################
#
# Update CanopyLayer struct when temperature or upstream pressure changes
#
###############################################################################
"""
    update_leaf_TP!(
                photo_set::AbstractPhotoModelParaSet{FT},
                canopyi::CanopyLayer{FT},
                hs::LeafHydraulics{FT},
                envir::AirLayer{FT}
    ) where {FT<:AbstractFloat}

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
        canopyi.g_max = g_max25 * relative_diffusive_coefficient(T);
        canopyi.g_min = g_min25 * relative_diffusive_coefficient(T);
        canopyi.LV    = latent_heat_vapor(T) * M_H₂O(FT);
        canopyi.ps.T  = canopyi.T;
        leaf_temperature_dependence!(photo_set, canopyi.ps, envir);
        canopyi.p_sat = canopyi.ps.p_sat;
        temperature_effects!(hs);
        hs.p_ups      = canopyi.p_ups;
        canopyi.ec    = critical_flow(hs, canopyi.ec);
        canopyi.T_old = canopyi.T;
        canopyi.p_old = canopyi.p_ups;
    # if only p_ups changes, update ec and p_old
    elseif p_ups != p_old
        hs.p_ups      = canopyi.p_ups;
        canopyi.ec    = critical_flow(hs, canopyi.ec);
        canopyi.p_old = canopyi.p_ups;
    end

    return nothing
end




function update_leaf_TP!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            hs::TreeSimple{FT},
            envir::AirLayer{FT}
) where {FT<:AbstractFloat}
    # unpack required variables
    @unpack g_max25, g_min25, p_ups, p_old, T, T_old = canopyi;

    # note that canopyi.p_ups is not used here
    # if T changes, update TD and ec, then T_old
    if T != T_old
        canopyi.g_max = g_max25 * relative_diffusive_coefficient(T);
        canopyi.g_min = g_min25 * relative_diffusive_coefficient(T);
        canopyi.LV    = latent_heat_vapor(T) * M_H₂O(FT);
        canopyi.ps.T  = canopyi.T;
        leaf_temperature_dependence!(photo_set, canopyi.ps, envir);
        canopyi.p_sat = canopyi.ps.p_sat;
        temperature_effects!(hs);
        tree_ec       = critical_flow(hs, canopyi.ec * canopyi.LA);
        canopyi.ec    = tree_ec / canopyi.LA;
        canopyi.T_old = canopyi.T;
        canopyi.p_old = canopyi.p_ups;
    # if only p_ups changes, update ec and p_old
    elseif p_ups != p_old
        tree_ec       = critical_flow(hs, canopyi.ec * canopyi.LA);
        canopyi.ec    = tree_ec / canopyi.LA;
        canopyi.p_old = canopyi.p_ups;
    end

    return nothing
end








###############################################################################
#
# Update leaf maximal A and K particularly for optimization models
#
###############################################################################
"""
    update_leaf_AK!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            envir::AirLayer{FT}
    ) where {FT<:AbstractFloat}

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
    for i in eachindex(APAR)
        _g_sw = 1 / max(1/_g_crit - 1/g_bw[i], FT(1e-3));
        _g_sw = min(_g_sw, g_max);
        _g_sw = max(_g_sw, g_min);
        _g_lc = 1 / (1/g_bc[i] + FT(1.6)/_g_sw + 1/g_m[i]);
        canopyi.ps.APAR = APAR[i];
        leaf_photosynthesis!(photo_set, canopyi.ps, envir, GCO₂Mode(), _g_lc);

        # update the a_max for each leaf
        canopyi.a_max[i] = canopyi.ps.An;
    end

    # update the kr_max
    canopyi.kr_max  = hs.k_history[end];

    return nothing
end
