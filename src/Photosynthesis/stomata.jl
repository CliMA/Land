###############################################################################
#
# Calculate empirical gsw from the equations
#
###############################################################################
"""
    empirical_gsw_from_model(
            model::EmpiricalStomatalModel,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            β::FT)

Steady state gsw from empirical approach given
- `model` [`EmpiricalStomatalModel`](@ref) type empirical model parameter set
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
- `β` Correction factor over the g1 part of an empirical model
"""
function empirical_gsw_from_model(
            model::ESMBallBerry{FT},
            leaf::AbstractLeaf{FT},
            envir::AirLayer{FT},
            β::FT
            ) where {FT<:AbstractFloat}
    @unpack g0, g1    = model;
    @unpack An, p_s   = leaf;
    @unpack p_atm, RH = envir;

    return g0 .+ g1 * RH * p_atm * FT(1e-6) * β .* An ./ p_s
end

function empirical_gsw_from_model(
            model::ESMGentine{FT},
            leaf::AbstractLeaf{FT},
            envir::AirLayer{FT},
            β::FT
            ) where {FT<:AbstractFloat}
    @unpack g0, g1  = model;
    @unpack An, p_i = leaf;
    @unpack p_atm   = envir;

    return g0 .+ g1 * p_atm * FT(1e-6) * β .* An ./ p_i
end

function empirical_gsw_from_model(
            model::ESMLeuning{FT},
            leaf::AbstractLeaf{FT},
            envir::AirLayer{FT},
            β::FT
            ) where {FT<:AbstractFloat}
    @unpack d0, g0, g1             = model;
    @unpack An, p_s, p_sat, Γ_star = leaf;
    @unpack p_atm, p_H₂O           = envir;

    return g0 .+ g1 * p_atm * FT(1e-6) / (1 + (p_sat - p_H₂O)/d0) *
                 β .* An ./ (p_s .- Γ_star)
                    
end

function empirical_gsw_from_model(
            model::ESMMedlyn{FT},
            leaf::AbstractLeaf{FT},
            envir::AirLayer{FT},
            β::FT
            ) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack An, p_sat         = leaf;
    @unpack p_a, p_atm, p_H₂O = envir;

    return g0 .+ p_atm * FT(1e-6) / p_a * (1 + g1/sqrt(p_sat - p_H₂O)) *
                 β .* An
end








###############################################################################
#
# Make g_min and g_max control
#
###############################################################################
"""
    leaf_gsw_control!(
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT})

make sure g_sw is in its physiological range limited by diffusion, given
- `photo_set` [`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type parameter set
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
"""
function leaf_gsw_control!(
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT}
            ) where {FT<:AbstractFloat}
    # if gsw is too small, use g_min
    if leaf.g_sw < leaf.g_min
        leaf.g_sw = leaf.g_min;
        leaf.g_sc = leaf.g_sw / FT(1.6);
        leaf.g_lw = 1 / (1/leaf.g_sw + 1/leaf.g_bw);
        leaf.g_lc = 1 / (1/leaf.g_bc + 1/leaf.g_sc + 1/leaf.g_m);

        leaf_photo_from_glc!(photo_set, leaf, envir);
    # if gsw is too big, use g_max
    elseif leaf.g_sw > leaf.g_max
        leaf.g_sw = leaf.g_max;
        leaf.g_sc = leaf.g_sw / FT(1.6);
        leaf.g_lw = 1 / (1/leaf.g_sw + 1/leaf.g_bw);
        leaf.g_lc = 1 / (1/leaf.g_bc + 1/leaf.g_sc + 1/leaf.g_m);

        leaf_photo_from_glc!(photo_set, leaf, envir);
    end

    # update leaf flow rate
    leaf.e = leaf.g_lw * (leaf.p_sat - envir.p_H₂O) / envir.p_atm;

    return nothing
end

function leaf_gsw_control!(
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaves{FT},
            envir::AirLayer{FT}
            ) where {FT<:AbstractFloat}
    count = 0
    for i in eachindex(leaf.g_sw)
        # if gsw is too small, use g_min
        if leaf.g_sw[i] < leaf.g_min
            leaf.g_sw[i] = leaf.g_min;
            leaf.g_sc[i] = leaf.g_sw[i] / FT(1.6);
            leaf.g_lw[i] = 1 / (1/leaf.g_sw[i] +
                                1/leaf.g_bw[i]);
            leaf.g_lc[i] = 1 / (1/leaf.g_bc[i] +
                                1/leaf.g_sc[i] +
                                1/leaf.g_m[i]);
            count += 1;
        # if gsw is too big, use g_max
        elseif leaf.g_sw[i] < leaf.g_min
            leaf.g_sw[i] = leaf.g_max;
            leaf.g_sc[i] = leaf.g_sw[i] / FT(1.6);
            leaf.g_lw[i] = 1 / (1/leaf.g_sw[i] +
                                1/leaf.g_bw[i]);
            leaf.g_lc[i] = 1 / (1/leaf.g_bc[i] +
                                1/leaf.g_sc[i] +
                                1/leaf.g_m[i]);
            count += 1;
        end
    end

    # update the photosynthetic rates if count > 0
    if count > 0
        leaf_photo_from_glc!(photo_set, leaf, envir);
    end

    # update leaf flow rate
    _d      = (leaf.p_sat - envir.p_H₂O) / envir.p_atm
    leaf.e .= leaf.g_lw * _d;

    return nothing
end








###############################################################################
#
# Diff function to minimize by RootSolvers
#
###############################################################################
"""
    envir_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::AbstractStomatalModel)

Calculate the difference to be minimized for a given
- `x` Assumed leaf diffusive conductance
- `photo_set`[`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type parameter set
- `leaf`[`Leaf`](@ref) type struct
- `envir`[`AirLayer`](@ref) type struct
- `sm` Stomatal model option (photo_set.Sto)
"""
function envir_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::ESMGentine{FT}
            ) where {FT<:AbstractFloat}
    # update photosynthesis from g_lc
    leaf.g_lc = x;
    leaf_photo_from_glc!(photo_set, leaf, envir);

    # unpack unchanged variables from structs
    @unpack An, g_bc, g_bw, g_lc, g_m, hs, p_sat = leaf;
    @unpack p_a, p_atm, p_H₂O = envir;

    # calculate flow rate
    g_sc = 1 / max( 1/g_lc - 1/g_m - 1/g_bc, FT(1.6e-3) );
    g_sw = g_sc * FT(1.6);
    g_lw = 1 / (1/g_sw + 1/g_bw);
    e_lf = g_lw * (p_sat - p_H₂O) / p_atm;

    # calculate the xylem end pressure
    k_lf = leaf_xylem_risk(hs, e_lf);

    # calculate g_sw from stomatal model
    g_mod = empirical_gsw_from_model(sm, leaf, envir, k_lf);
    g_mod = min(leaf.g_max, g_mod);

    # calculate model predicted p_i
    g_leaf = 1 / (FT(1.6)/g_mod + 1/leaf.g_bc + 1/leaf.g_m);

    return g_leaf - x
end

function envir_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaves{FT},
            envir::AirLayer{FT},
            sm::ESMGentine{FT},
            ind::Int
            ) where {FT<:AbstractFloat}
    # update photosynthesis from g_lc
    leaf.g_lc[ind] = x;
    leaf_photo_from_glc!(photo_set, leaf, envir, ind);

    # unpack unchanged variables from structs
    @unpack An, g_bc, g_bw, g_lc, g_m, hs, p_sat = leaf;
    @unpack p_a, p_atm, p_H₂O = envir;

    # calculate flow rate
    g_sc = 1 / max( 1/x - 1/g_m[ind] - 1/g_bc[ind], FT(1.6e-3) );
    g_sw = g_sc * FT(1.6);
    g_lw = 1 / (1/g_sw + 1/g_bw);
    e_lf = g_lw * (p_sat - p_H₂O) / p_atm;

    # calculate the xylem end pressure
    k_lf = leaf_xylem_risk(hs, e_lf);

    # calculate g_sw from stomatal model
    g_mod = empirical_gsw_from_model(sm, leaf, envir, k_lf);
    g_mod = min(leaf.g_max, g_mod);

    # calculate model predicted p_i
    g_leaf = 1 / (FT(1.6)/g_mod + 1/g_bc[ind] + 1/g_m[ind]);

    return g_leaf - x
end

function envir_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::EmpiricalStomatalModel
            ) where {FT<:AbstractFloat}
    # update photosynthesis from g_lc
    leaf.g_lc = x;
    leaf_photo_from_glc!(photo_set, leaf, envir);

    # calculate g_sw from stomatal model
    g_mod = empirical_gsw_from_model(sm, leaf, envir, FT(1));
    g_mod = min(leaf.g_max, g_mod);

    # calculate model predicted p_i
    g_leaf = 1 / (FT(1.6)/g_mod + 1/leaf.g_bc + 1/leaf.g_m);

    return g_leaf - x
end

function envir_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaves{FT},
            envir::AirLayer{FT},
            sm::EmpiricalStomatalModel,
            ind::Int
            ) where {FT<:AbstractFloat}
    # update photosynthesis from g_lc
    leaf.g_lc[ind] = x;
    leaf_photo_from_glc!(photo_set, leaf, envir, ind);

    # calculate g_sw from stomatal model
    g_mod = empirical_gsw_from_model(sm, leaf, envir, FT(1));
    g_mod = min(leaf.g_max, g_mod);

    # calculate model predicted p_i
    g_leaf = 1 / (FT(1.6)/g_mod + 1/leaf.g_bc[ind] + 1/leaf.g_m[ind]);

    return g_leaf - x
end

function envir_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::OSMEller
            ) where {FT<:AbstractFloat}
    # calculate the limit of g_lc to ensure x in the physiological range
    _glh = 1 / (1/leaf.g_bc + FT(1.6)/leaf.g_max + 1/leaf.g_m);
    _gll = 1 / (1/leaf.g_bc + FT(1.6)/leaf.g_min + 1/leaf.g_m);

    if x > _glh
        x = _glh + FT(1e-4);

        return FT(0)
    elseif x< _gll
        x = _gll - FT(1e-4);

        return FT(0)
    else
        # update photosynthesis from g_lc-FT(1e-3)
        leaf.g_lc = x-FT(1e-3);
        leaf_photo_from_glc!(photo_set, leaf, envir);

        @unpack An, g_bc, g_bw, g_lc, g_m, hs, p_sat = leaf;
        @unpack p_atm, p_H₂O = envir;

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        a_1  = An;
        e_1  = g_lw * (p_sat - p_H₂O) / p_atm;
        k_1  = leaf_xylem_risk(hs, e_1);

        # update photosynthesis from g_lc
        leaf.g_lc = x;
        leaf_photo_from_glc!(photo_set, leaf, envir);

        @unpack An, g_bc, g_bw, g_lc, g_m, hs, p_sat = leaf;
        @unpack p_atm, p_H₂O = envir;
        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        a_2  = An;
        e_2  = g_lw * (p_sat - p_H₂O) / p_atm;
        k_2  = leaf_xylem_risk(hs, e_2);

        ∂A∂E = (a_2 - a_1) / (e_2 - e_1);
        ∂Θ∂E = (k_1 - k_2) / (e_2 - e_1) * a_2 / k_2;
        diff = ∂A∂E - ∂Θ∂E;

        #= used for debugging
        @show x;
        @show ∂A∂E;
        @show ∂Θ∂E;
        println("");
        #sleep(0.1);
        # =#

        return diff
    end
end

function envir_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::OSMSperry
            ) where {FT<:AbstractFloat}
    # calculate the limit of g_lc to ensure x in the physiological range
    _glh = 1 / (1/leaf.g_bc + FT(1.6)/leaf.g_max + 1/leaf.g_m);
    _gll = 1 / (1/leaf.g_bc + FT(1.6)/leaf.g_min + 1/leaf.g_m);

    if x > _glh
        x = _glh + FT(1e-4);

        return FT(0)
    elseif x< _gll
        x = _gll - FT(1e-4);

        return FT(0)
    else
        # update photosynthesis from g_lc-FT(1e-3)
        leaf.g_lc = x-FT(1e-3);
        leaf_photo_from_glc!(photo_set, leaf, envir);

        @unpack a_max, An, g_bc, g_bw, g_lc, g_m, hs, kr_max, p_sat = leaf;
        @unpack p_atm, p_H₂O = envir;

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        a_1  = An;
        e_1  = g_lw * (p_sat - p_H₂O) / p_atm;
        k_1  = leaf_xylem_risk(hs, e_1);

        # update photosynthesis from g_lc
        leaf.g_lc = x;
        leaf_photo_from_glc!(photo_set, leaf, envir);

        @unpack a_max, An, g_bc, g_bw, g_lc, g_m, hs, kr_max, p_sat = leaf;
        @unpack p_atm, p_H₂O = envir;

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        a_2  = An;
        e_2  = g_lw * (p_sat - p_H₂O) / p_atm;
        k_2  = leaf_xylem_risk(hs, e_2);

        ∂A∂E = (a_2 - a_1) / (e_2 - e_1);
        ∂Θ∂E = (k_1 - k_2) / (e_2 - e_1) * a_max / kr_max;
        diff = ∂A∂E - ∂Θ∂E;

        #= used for debugging
        @show x;
        @show ∂A∂E;
        @show ∂Θ∂E;
        println("");
        #sleep(0.1);
        # =#

        return diff
    end
end

function envir_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::OSMWang
            ) where {FT<:AbstractFloat}
    # calculate the limit of g_lc to ensure x in the physiological range
    _glh = 1 / (1/leaf.g_bc + FT(1.6)/leaf.g_max + 1/leaf.g_m);
    _gll = 1 / (1/leaf.g_bc + FT(1.6)/leaf.g_min + 1/leaf.g_m);

    if x > _glh
        x = _glh + FT(1e-4);

        return FT(0)
    elseif x< _gll
        x = _gll - FT(1e-4);

        return FT(0)
    else
        # update photosynthesis from g_lc-FT(1e-3)
        leaf.g_lc = x-FT(1e-3);
        leaf_photo_from_glc!(photo_set, leaf, envir);

        @unpack An, ec, g_bc, g_bw, g_lc, g_m, p_sat = leaf;
        @unpack p_atm, p_H₂O = envir;

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        e_1  = g_lw * (p_sat - p_H₂O) / p_atm;
        a_1  = An;

        # update photosynthesis from g_lc
        leaf.g_lc = x;
        leaf_photo_from_glc!(photo_set, leaf, envir);

        @unpack An, ec, g_bc, g_bw, g_lc, g_m, p_sat = leaf;
        @unpack p_atm, p_H₂O = envir;

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        e_2  = g_lw * (p_sat - p_H₂O) / p_atm;
        a_2  = An;

        ∂A∂E = (a_2 - a_1) / (e_2 - e_1);
        ∂Θ∂E = a_2 / max(ec - e_2, FT(1e-7));
        diff = ∂A∂E - ∂Θ∂E;

        #= used for debugging
        @show x;
        @show ∂A∂E;
        @show ∂Θ∂E;
        println("");
        #sleep(0.1);
        # =#

        return diff
    end
end

function envir_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::OSMWAP
            ) where {FT<:AbstractFloat}
    # calculate the limit of g_lc to ensure x in the physiological range
    _glh = 1 / (1/leaf.g_bc + FT(1.6)/leaf.g_max + 1/leaf.g_m);
    _gll = 1 / (1/leaf.g_bc + FT(1.6)/leaf.g_min + 1/leaf.g_m);

    if x > _glh
        x = _glh + FT(1e-4);

        return FT(0)
    elseif x< _gll
        x = _gll - FT(1e-4);

        return FT(0)
    else
        # update photosynthesis from g_lc-FT(1e-3)
        leaf.g_lc = x-FT(1e-3);
        leaf_photo_from_glc!(photo_set, leaf, envir);

        @unpack An, g_bc, g_bw, g_lc, g_m, hs, p_sat = leaf;
        @unpack p_a, p_atm, p_H₂O = envir;

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        a_1  = An;
        e_1  = g_lw * (p_sat - p_H₂O) / p_atm;
        p_1  = xylem_p_from_flow(hs, e_1);

        # update photosynthesis from g_lc-FT(1e-4)
        leaf.g_lc = x;
        leaf_photo_from_glc!(photo_set, leaf, envir);

        @unpack An, g_bc, g_bw, g_lc, g_m, hs, p_sat = leaf;
        @unpack p_a, p_atm, p_H₂O = envir;

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        a_2  = An;
        e_2  = g_lw * (p_sat - p_H₂O) / p_atm;
        p_2  = xylem_p_from_flow(hs, e_2);

        ∂A∂E = (a_2 - a_1) / (e_2 - e_1);
        ∂Θ∂E = (-2 * sm.a * p_2 + sm.b) / (e_2 - e_1) * (p_1 - p_2);
        diff = ∂A∂E - ∂Θ∂E;

        #= used for debugging
        @show x;
        @show ∂A∂E;
        @show ∂Θ∂E;
        println("");
        #sleep(0.1);
        # =#

        return diff
    end
end

function envir_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::OSMWAPMod
            ) where {FT<:AbstractFloat}
    # calculate the limit of g_lc to ensure x in the physiological range
    _glh = 1 / (1/leaf.g_bc + FT(1.6)/leaf.g_max + 1/leaf.g_m);
    _gll = 1 / (1/leaf.g_bc + FT(1.6)/leaf.g_min + 1/leaf.g_m);

    if x > _glh
        x = _glh + FT(1e-4);

        return FT(0)
    elseif x< _gll
        x = _gll - FT(1e-4);

        return FT(0)
    else
        # update photosynthesis from g_lc-FT(1e-3)
        leaf.g_lc = x-FT(1e-3);
        leaf_photo_from_glc!(photo_set, leaf, envir);

        @unpack An, g_bc, g_bw, g_lc, g_m, hs, p_sat = leaf;
        @unpack p_a, p_atm, p_H₂O = envir;

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        a_1  = An;
        e_1  = g_lw * (p_sat - p_H₂O) / p_atm;
        p_1  = xylem_p_from_flow(hs, e_1);

        # update photosynthesis from g_lc-FT(1e-4)
        leaf.g_lc = x;
        leaf_photo_from_glc!(photo_set, leaf, envir);

        @unpack An, g_bc, g_bw, g_lc, g_m, hs, p_sat = leaf;
        @unpack p_a, p_atm, p_H₂O = envir;

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        a_2  = An;
        e_2  = g_lw * (p_sat - p_H₂O) / p_atm;
        p_2  = xylem_p_from_flow(hs, e_2);

        ∂A∂E = (a_2 - a_1) / (e_2 - e_1);
        ∂Θ∂E = (-2 * sm.a * p_2 * a_2) / (e_2 - e_1) * (p_1 - p_2);
        diff = ∂A∂E - ∂Θ∂E;

        #= used for debugging
        @show x;
        @show ∂A∂E;
        @show ∂Θ∂E;
        println("");
        #sleep(0.1);
        # =#

        return diff
    end
end








###############################################################################
#
# Solve for stomatal conductance from environmental conditions
#
###############################################################################
"""
    leaf_photo_from_envir!(
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::AbstractStomatalModel)

Calculate steady state gsw and photosynthesis from empirical approach, given
- `photo_set` [`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type parameter set
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
- `sm` [`EmpiricalStomatalModel`](@ref) or [`OptimizationStomatalModel`](@ref)

Note that `sm` is stored in `photo_set.Sto`
"""
function leaf_photo_from_envir!(
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::OSMEller
            ) where {FT<:AbstractFloat}
    if leaf.APAR > 1
        _hs::LeafHydraulics = leaf.hs

        # if T changes, update TD and ec, then T_old
        if leaf.T != leaf.T_old
            leaf_temperature_dependence!(photo_set, leaf, envir);
            leaf.ec = leaf_e_crit(_hs);
            leaf.T_old = leaf.T;
        # if only p_ups changes, update ec and p_old
        elseif leaf.p_ups != leaf.p_old
            leaf.ec = leaf_e_crit(_hs);
            leaf.p_old = leaf.p_ups;
        end

        # calculate the physiological maximal g_sw
        # add an 90% offset in _g_crit for Eller model to avoid dK/dE = 0
        _g_crit = FT(0.9) * leaf.ec / (leaf.p_sat - envir.p_H₂O) * envir.p_atm;
        _g_sw   = 1 / max(1/_g_crit - 1/leaf.g_bw, FT(1e-3));

        # if _g_sw is lower than g_min
        if _g_sw <= leaf.g_min
            leaf.g_sw = 0

        # if _g_sw is higher than g_min, reculate A_max and Kr_max
        else
            _g_max = min(_g_sw, leaf.g_max);
            leaf.g_sw   = _g_max;
            leaf.g_lw   = 1 / (1/leaf.g_bw + 1/leaf.g_sw);
            leaf.g_sc   = leaf.g_sw / FT(1.6);
            leaf.g_lc   = 1 / (1/leaf.g_bc + 1/leaf.g_sc + 1/leaf.g_m);
            leaf_photo_from_glc!(photo_set, leaf, envir);
            leaf.a_max  = leaf.An;
            leaf.kr_max = _hs.k_history[end];

            # if net A is too small
            if leaf.a_max < 0.05
                leaf.g_sw = FT(0)

            # if net A is high enough
            else
                # solve for optimal g_lc, A and g_sw updated here
                _gh    = 1 / (1/leaf.g_bc + FT(1.6)/_g_max + 1/leaf.g_m);
                _gl    = 1 / (1/leaf.g_bc + FT(1.6)/leaf.g_min + 1/leaf.g_m);
                _sm    = NewtonBisectionMethod{FT}(_gl, _gh);
                _st    = SolutionTolerance{FT}(1e-4);
                @inline f(x) = envir_diff!(x, photo_set, leaf, envir, sm);
                _solut = find_zero_ext(f, _sm, _st);

                #= used for debugging
                @show _g_sw;
                @show _g_max;
                @show leaf.p_ups;
                @show _solut;
                println("");
                =#

                # update leaf conductances
                leaf.g_lc = _solut;
                leaf.g_sc = 1 / (1/leaf.g_lc - 1/leaf.g_m - 1/leaf.g_bc);
                leaf.g_sw = leaf.g_sc * FT(1.6);
                leaf.g_lw = 1 / (1/leaf.g_sw + 1/leaf.g_bw);
                leaf_photo_from_glc!(photo_set, leaf, envir);
            end
        end
    else
        leaf.g_sw = 0;
    end

    # make sure g_sw in its range, and update flow rate
    leaf_gsw_control!(photo_set, leaf, envir);

    return nothing
end

function leaf_photo_from_envir!(
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::OSMSperry
            ) where {FT<:AbstractFloat}
    if leaf.APAR > 1
        _hs::LeafHydraulics = leaf.hs

        # if T changes, update TD and ec, then T_old
        if leaf.T != leaf.T_old
            leaf_temperature_dependence!(photo_set, leaf, envir);
            leaf.ec = leaf_e_crit(_hs);
            leaf.T_old = leaf.T;
        # if only p_ups changes, update ec and p_old
        elseif leaf.p_ups != leaf.p_old
            leaf.ec = leaf_e_crit(_hs);
            leaf.p_old = leaf.p_ups;
        end

        # calculate the physiological maximal g_sw
        # add an 70% offset in _g_crit for Sperry model to avoid dK/dE = 0
        _g_crit = FT(0.7) * leaf.ec / (leaf.p_sat - envir.p_H₂O) * envir.p_atm;
        _g_sw   = 1 / max(1/_g_crit - 1/leaf.g_bw, FT(1e-3));

        # if _g_sw is lower than g_min
        if _g_sw <= leaf.g_min
            leaf.g_sw = 0

        # if _g_sw is higher than g_min, reculate A_max and Kr_max
        else
            _g_max = min(_g_sw, leaf.g_max);
            leaf.g_sw   = _g_max;
            leaf.g_lw   = 1 / (1/leaf.g_bw + 1/leaf.g_sw);
            leaf.g_sc   = leaf.g_sw / FT(1.6);
            leaf.g_lc   = 1 / (1/leaf.g_bc + 1/leaf.g_sc + 1/leaf.g_m);
            leaf_photo_from_glc!(photo_set, leaf, envir);
            leaf.a_max  = leaf.An;
            leaf.kr_max = _hs.k_history[end];

            # if net A is too small
            if leaf.a_max < 0.05
                leaf.g_sw = FT(0)

            # if net A is high enough
            else
                # solve for optimal g_lc, A and g_sw updated here
                _gh    = 1 / (1/leaf.g_bc + FT(1.6)/_g_max + 1/leaf.g_m);
                _gl    = 1 / (1/leaf.g_bc + FT(1.6)/leaf.g_min + 1/leaf.g_m);
                _sm    = NewtonBisectionMethod{FT}(_gl, _gh);
                _st    = SolutionTolerance{FT}(1e-4);
                @inline f(x) = envir_diff!(x, photo_set, leaf, envir, sm);
                _solut = find_zero_ext(f, _sm, _st);

                #= used for debugging
                @show _g_sw;
                @show _g_max;
                @show leaf.p_ups;
                @show _solut;
                println("");
                =#

                # update leaf conductances
                leaf.g_lc = _solut;
                leaf.g_sc = 1 / (1/leaf.g_lc - 1/leaf.g_m - 1/leaf.g_bc);
                leaf.g_sw = leaf.g_sc * FT(1.6);
                leaf.g_lw = 1 / (1/leaf.g_sw + 1/leaf.g_bw);
                leaf_photo_from_glc!(photo_set, leaf, envir);
            end
        end
    else
        leaf.g_sw = 0;
    end

    # make sure g_sw in its range, and update flow rate
    leaf_gsw_control!(photo_set, leaf, envir);

    return nothing
end

function leaf_photo_from_envir!(
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::OptimizationStomatalModel
            ) where {FT<:AbstractFloat}
    if leaf.APAR > 1
        _hs::LeafHydraulics = leaf.hs

        # if T changes, update TD and ec, then T_old
        if leaf.T != leaf.T_old
            leaf_temperature_dependence!(photo_set, leaf, envir);
            leaf.ec = leaf_e_crit(_hs);
            leaf.T_old = leaf.T;
        # if only p_ups changes, update ec and p_old
        elseif leaf.p_ups != leaf.p_old
            leaf.ec = leaf_e_crit(_hs);
            leaf.p_old = leaf.p_ups;
        end

        # calculate the physiological maximal g_sw
        _g_crit = leaf.ec / (leaf.p_sat - envir.p_H₂O) * envir.p_atm;
        _g_sw   = 1 / max(1/_g_crit - 1/leaf.g_bw, FT(1e-3));

        # if _g_sw is lower than g_min
        if _g_sw <= leaf.g_min
            leaf.g_sw = 0

        # if _g_sw is higher than g_min, reculate A_max and Kr_max
        else
            _g_max = min(_g_sw, leaf.g_max);
            leaf.g_sw   = _g_max;
            leaf.g_lw   = 1 / (1/leaf.g_bw + 1/leaf.g_sw);
            leaf.g_sc   = leaf.g_sw / FT(1.6);
            leaf.g_lc   = 1 / (1/leaf.g_bc + 1/leaf.g_sc + 1/leaf.g_m);
            leaf_photo_from_glc!(photo_set, leaf, envir);
            leaf.a_max  = leaf.An;
            leaf.kr_max = _hs.k_history[end];

            # if net A is too small
            if leaf.a_max < 0.05
                leaf.g_sw = FT(0)

            # if net A is high enough
            else
                # solve for optimal g_lc, A and g_sw updated here
                _gh    = 1 / (1/leaf.g_bc + FT(1.6)/_g_max + 1/leaf.g_m);
                _gl    = 1 / (1/leaf.g_bc + FT(1.6)/leaf.g_min + 1/leaf.g_m);
                _sm    = NewtonBisectionMethod{FT}(_gl, _gh);
                _st    = SolutionTolerance{FT}(1e-4);
                @inline f(x) = envir_diff!(x, photo_set, leaf, envir, sm);
                _solut = find_zero_ext(f, _sm, _st);

                #= used for debugging
                @show _g_sw;
                @show _g_max;
                @show leaf.p_ups;
                @show _solut;
                println("");
                =#

                # update leaf conductances
                leaf.g_lc = _solut;
                leaf.g_sc = 1 / (1/leaf.g_lc - 1/leaf.g_m - 1/leaf.g_bc);
                leaf.g_sw = leaf.g_sc * FT(1.6);
                leaf.g_lw = 1 / (1/leaf.g_sw + 1/leaf.g_bw);
                leaf_photo_from_glc!(photo_set, leaf, envir);
            end
        end
    else
        leaf.g_sw = 0;
    end

    # make sure g_sw in its range, and update flow rate
    leaf_gsw_control!(photo_set, leaf, envir);

    return nothing
end

function leaf_photo_from_envir!(
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::EmpiricalStomatalModel
            ) where {FT<:AbstractFloat}
    if leaf.APAR > 1
        # solve for optimal g_lc, A and g_sw updated here
        _gh = 1 / (1/leaf.g_bc + FT(1.6)/leaf.g_max + 1/leaf.g_m) - FT(1e-6);
        _gl = _gh - FT(0.001);
        _sm = SecantMethod{FT}(_gl, _gh);
        _cs = CompactSolution();
        _st = SolutionTolerance{FT}(1e-7);
        @inline f(x) = envir_diff!(x, photo_set, leaf, envir, sm);
        _solut = find_zero(f, _sm, _cs, _st);

        # update leaf conductances
        leaf.g_sc = 1 / max(1/leaf.g_lc - 1/leaf.g_m - 1/leaf.g_bc, FT(1.6e-3));
        leaf.g_sw = leaf.g_sc * FT(1.6);
        leaf.g_lw = 1 / (1/leaf.g_sw + 1/leaf.g_bw);
    else
        leaf.g_sw = 0;
    end

    # make sure g_sw in its range, and update flow rate
    leaf_gsw_control!(photo_set, leaf, envir);

    return nothing
end
