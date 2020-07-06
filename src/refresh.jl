###############################################################################
#
# Update Leaves struct when temperature changes
#
###############################################################################
"""
    update_leaf!(
            photo_set::AbstractPhotoModelParaSet{FT},
            leaves::Leaves{FT},
            envir::AirLayer{FT})

Update leaf physiological parameters if temperature or pressure changes in the
daytime, given
- `photo_set` [`C3ParaSet`] or [`C4ParaSet`] type parameter set
- `leaves` [`Leaves`](@ref) type struct
- `envir` [`AirLayer`] type struct
"""
function update_leaf_TP!(
            photo_set::AbstractPhotoModelParaSet{FT},
            leaves::Leaves{FT},
            envir::AirLayer{FT}
            ) where {FT<:AbstractFloat}
    # unpack required variables
    @unpack g_max25, g_min25, p_ups, p_old, T, T_old = leaves;

    # if T changes, update TD and ec, then T_old
    if T != T_old
        leaves.g_max    = g_max25 * relative_diffusive_coefficient(T);
        leaves.g_min    = g_min25 * relative_diffusive_coefficient(T);
        leaves.LV       = latent_heat_vapor(T) * 1000 / FT(MOLMASS_WATER);
        leaves.ps.T     = leaves.T;
        leaf_temperature_dependence!(photo_set, leaves.ps, envir);
        leaves.p_sat    = leaves.ps.p_sat;
        leaves.hs.f_st  = relative_surface_tension(T);
        leaves.hs.f_vis = relative_viscosity(T);
        leaves.hs.p_ups = leaves.p_ups;
        leaves.ec       = leaf_e_crit(leaves.hs, leaves.ec);
        leaves.T_old    = leaves.T;
        leaves.p_old    = leaves.p_ups;
    # if only p_ups changes, update ec and p_old
    elseif p_ups != p_old
        leaves.hs.p_ups = leaves.p_ups;
        leaves.ec       = leaf_e_crit(leaves.hs, leaves.ec);
        leaves.p_old    = leaves.p_ups;
    end

    return nothing
end




"""
    update_leaf_AK!(
            photo_set::AbstractPhotoModelParaSet{FT},
            leaves::Leaves{FT},
            envir::AirLayer{FT})

Update leaf maximal A and K for Sperry model, given
- `photo_set` [`C3ParaSet`] or [`C4ParaSet`] type parameter set
- `leaves` [`Leaves`](@ref) type struct
- `envir` [`AirLayer`] type struct
"""
function update_leaf_AK!(
            photo_set::AbstractPhotoModelParaSet{FT},
            leaves::Leaves{FT},
            envir::AirLayer{FT}
            ) where {FT<:AbstractFloat}
    # unpack required variables
    @unpack APAR, ec, g_bc, g_bw, g_m, g_max, g_min, p_sat = leaves;
    @unpack p_atm, p_H₂O = envir;

    # calculate the physiological maximal g_sw
    _g_crit = ec / (p_sat - p_H₂O) * p_atm;
    _g_sw   = 1 ./ max.(1/_g_crit .- 1 ./ g_bw, FT(1e-3));
    _g_sw  .= min.(_g_sw, g_max);
    _g_sw  .= max.(_g_sw, g_min);
    _g_lcs  = 1 ./ (1 ./ g_bc .+ FT(1.6) ./ _g_sw .+ 1 ./ g_m);

    # update the limited rates
    _Jps = leaf_ETR_pot_APAR(leaves.ps, APAR);
    _Js  = leaf_ETR_Jps(photo_set, leaves.ps, _Jps);
    _Anc = rubisco_limited_an_glc(photo_set, leaves.ps, envir, _g_lcs);
    _Anj = light_limited_an_glc(photo_set, leaves.ps, envir, _g_lcs, _Js);
    _Anp = product_limited_an_glc(photo_set, leaves.ps, envir, _g_lcs);
    _Ams = min.(_Anc, _Anj, _Anp);

    # update a_max and kr_max
    leaves.a_max  .= _Ams;
    leaves.kr_max  = leaves.hs.k_history[end];

    return nothing
end








###############################################################################
#
# Update Leaves struct from stomatal conductance
#
###############################################################################
"""
    update_leaf_from_glc!(
            photo_set::AbstractPhotoModelParaSet{FT},
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            ind::Int,
            glc::FT)

Update Nth leaf photosynthesis, given
- `photo_set` [`C3ParaSet`] or [`C4ParaSet`] type parameter set
- `leaves` [`Leaves`](@ref) type struct
- `envir` [`AirLayer`] type struct
- `ind` Nth leaf
- `glc` Given leaf diffusive conductance
"""
function update_leaf_from_glc!(
            photo_set::AbstractPhotoModelParaSet{FT},
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            ind::Int,
            glc::FT
            ) where {FT<:AbstractFloat}
    # update the conductances
    leaves.g_lc[ind] = glc;
    leaves.g_sc[ind] = 1 / ( 1 / glc -
                             1 / leaves.g_m[ind] -
                             1 / leaves.g_bc[ind] );
    leaves.g_sw[ind] = leaves.g_sc[ind] * FT(1.6);
    leaves.g_lw[ind] = 1 / ( 1 / leaves.g_sw[ind] +
                             1 / leaves.g_bw[ind] );

    # update the photosynthetic rates
    if glc != leaves.ps.g_lc
        leaf_photo_from_glc!(photo_set, leaves.ps, envir, glc);
    end
    leaves.Ac[ind] = leaves.ps.Ac;
    leaves.Aj[ind] = leaves.ps.Aj;
    leaves.Ap[ind] = leaves.ps.Ap;
    leaves.Ag[ind] = leaves.ps.Ag;
    leaves.An[ind] = leaves.ps.An;

    # update the pressures
    leaves.p_i[ind] = leaves.ps.p_i;
    leaves.p_s[ind] = leaves.ps.p_s;

    return nothing
end




"""
    update_leaf_from_gsw!(
            photo_set::AbstractPhotoModelParaSet{FT},
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            ind::Int,
            glc::FT)

Update Nth leaf photosynthesis, given
- `photo_set` [`C3ParaSet`] or [`C4ParaSet`] type parameter set
- `leaves` [`Leaves`](@ref) type struct
- `envir` [`AirLayer`] type struct
- `ind` Nth leaf
- `gsw` Given stomatal conductance to H₂O
"""
function update_leaf_from_gsw!(
            photo_set::AbstractPhotoModelParaSet{FT},
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            ind::Int,
            gsw::FT
            ) where {FT<:AbstractFloat}
    # update the conductances
    leaves.g_sw[ind] = gsw;
    leaves.g_lw[ind] = 1 / ( 1 / gsw +
                             1 / leaves.g_bw[ind] );
    leaves.g_sc[ind] = gsw / FT(1.6);
    leaves.g_lc[ind] = 1 / ( FT(1.6) / gsw +
                             1 / leaves.g_m[ind] +
                             1 / leaves.g_bc[ind] );

    # update the photosynthetic rates
    leaf_photo_from_glc!(photo_set, leaves.ps, envir, leaves.g_lc[ind]);
    leaves.Ac[ind] = leaves.ps.Ac;
    leaves.Aj[ind] = leaves.ps.Aj;
    leaves.Ap[ind] = leaves.ps.Ap;
    leaves.Ag[ind] = leaves.ps.Ag;
    leaves.An[ind] = leaves.ps.An;

    # update the pressures
    leaves.p_i[ind] = leaves.ps.p_i;
    leaves.p_s[ind] = leaves.ps.p_s;

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
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            ind::Int)

make sure g_sw is in its physiological range limited by diffusion, given
- `photo_set` [`C3ParaSet`] or [`C4ParaSet`] type parameter set
- `leaves` [`Leaves`](@ref) type struct
- `envir` [`AirLayer`] type struct
"""
function leaf_gsw_control!(
            photo_set::AbstractPhotoModelParaSet,
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            ind::Int
            ) where {FT<:AbstractFloat}
    # if g_sw is low than g_min
    if leaves.g_sw[ind] < leaves.g_min
        update_leaf_from_gsw!(photo_set, leaves, envir, ind, leaves.g_min);
    elseif leaves.g_sw[ind] > leaves.g_max
        update_leaf_from_gsw!(photo_set, leaves, envir, ind, leaves.g_max);
    end

    return nothing
end
