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
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            β::FT
            ) where {FT<:AbstractFloat}
    @unpack g0, g1    = model;
    @unpack An, p_s   = leaf;
    @unpack p_atm, RH = envir;

    return g0 + β * g1 * RH * An / (p_s / p_atm * FT(1e6))
end

function empirical_gsw_from_model(
            model::ESMGentine{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            β::FT
            ) where {FT<:AbstractFloat}
    @unpack g0, g1  = model;
    @unpack An, p_i = leaf;
    @unpack p_atm   = envir;

    return g0 + β * g1 * An / (p_i / p_atm * FT(1e6))
end

function empirical_gsw_from_model(
            model::ESMLeuning{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            β::FT
            ) where {FT<:AbstractFloat}
    @unpack d0, g0, g1             = model;
    @unpack An, p_s, p_sat, Γ_star = leaf;
    @unpack p_atm, p_H₂O           = envir;

    return g0 + β * g1 * An / ((p_s - Γ_star) / p_atm * FT(1e6)) /
                (1 + (p_sat - p_H₂O)/d0)
end

function empirical_gsw_from_model(
            model::ESMMedlyn{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            β::FT
            ) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack An, p_sat         = leaf;
    @unpack p_a, p_atm, p_H₂O = envir;

    return g0 + β * (1 + g1/sqrt(p_sat - p_H₂O)) * An / (p_a / p_atm * FT(1e6))
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
    end

    # if gsw is too big, use g_max
    if leaf.g_sw > leaf.g_max
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
- `x` Assumed leaf internal CO₂
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
    # update photosynthesis from p_i
    leaf.p_i = x;
    photo_CO₂_dependence!(photo_set, leaf);

    # unpack unchanged variables from structs
    @unpack An, g_bc, g_bw, g_m, hs, p_sat = leaf;
    @unpack p_a, p_atm, p_H₂O = envir;

    # calculate flow rate
    g_lc = An*FT(1e-6) / (p_a - x) * p_atm
    g_sc = 1 / max( 1/g_lc - 1/g_m - 1/g_bc, FT(1e-3) )
    g_sw = g_sc * FT(1.6)
    g_lw = 1 / (1/g_sw + 1/g_bw)
    e_lf = g_lw * (p_sat - p_H₂O) / p_atm

    # calculate the xylem end pressure
    k_lf = leaf_xylem_risk(hs, e_lf);

    # calculate g_sw from stomatal model
    if leaf.An > 0
        g_mod = empirical_gsw_from_model(sm, leaf, envir, k_lf);
        g_mod = min(leaf.g_max, g_mod);
    else
        g_mod = leaf.g_min;
    end
    leaf.g_sw = g_mod;

    # calculate model predicted p_i
    g_leaf = 1 / (FT(1.6)/g_mod + 1/leaf.g_bc + 1/leaf.g_m);
    mod_x  = p_a - leaf.An*FT(1e-6) * p_atm / g_leaf;

    return mod_x - x
end

function envir_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::EmpiricalStomatalModel
            ) where {FT<:AbstractFloat}
    # update photosynthesis from p_i
    leaf.p_i = x;
    photo_CO₂_dependence!(photo_set, leaf);
    leaf.p_s = envir.p_a - leaf.An*FT(1e-6) * envir.p_atm / leaf.g_bc;

    # calculate g_sw from stomatal model
    if leaf.An > 0
        g_mod = empirical_gsw_from_model(sm, leaf, envir, FT(1));
        g_mod = min(leaf.g_max, g_mod);
    else
        g_mod = leaf.g_min;
    end

    # calculate model predicted p_i
    g_leaf = 1 / (FT(1.6)/g_mod + 1/leaf.g_bc + 1/leaf.g_m);
    mod_x  = envir.p_a - leaf.An*FT(1e-6) * envir.p_atm / g_leaf;

    return mod_x - x
end

function envir_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::OSMEller
            ) where {FT<:AbstractFloat}
    @unpack g_bc, g_bw, g_m, hs, p_sat = leaf;
    @unpack p_a, p_atm, p_H₂O = envir;

    # Eller model run into numerical issues when K = 1 using Float32
    # Thus, a control of maximal p_i = p_a applies
    if (x >= p_a) || (isnan(x))
        x = p_a;
        leaf.p_i = x;
        photo_CO₂_dependence!(photo_set, leaf);

        return FT(0)
    else
        # calculate flow rate and optimizer at the given p_i - Δp
        leaf.p_i = x - FT(1e-4);
        photo_CO₂_dependence!(photo_set, leaf);
        g_lc = leaf.An*FT(1e-6) / (p_a - leaf.p_i) * p_atm;
        g_sc = 1 / max( 1/g_lc - 1/g_m - 1/g_bc, FT(1e-3) );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        e_1  = g_lw * (p_sat - p_H₂O) / p_atm;
        k_1  = leaf_xylem_risk(hs, e_1);
        opt1 = leaf.An * k_1;
    
        # calculate flow rate and optimizer at the given p_i
        leaf.p_i = x;
        photo_CO₂_dependence!(photo_set, leaf);
        g_lc = leaf.An*FT(1e-6) / (p_a - leaf.p_i) * p_atm;
        g_sc = 1 / max( 1/g_lc - 1/g_m - 1/g_bc, FT(1e-3) );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        e_2  = g_lw * (p_sat - p_H₂O) / p_atm;
        k_2  = leaf_xylem_risk(hs, e_2);
        opt2 = leaf.An * k_2;
    
        # calculate the marginal gain - marginal penalty
        slope = (opt2 - opt1) / (e_2 - e_1) * FT(1e-3);
    
        return slope
    end
end

function envir_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::OSMSperry
            ) where {FT<:AbstractFloat}
    @unpack a_max, g_bc, g_bw, g_m, hs, kr_max, p_sat = leaf;
    @unpack p_a, p_atm, p_H₂O = envir;

    # Sperry model run into numerical issues when K = 1 using Float32
    # Thus, a control of maximal p_i = p_a applies
    if (x >= p_a) || (isnan(x))
        x = p_a;
        leaf.p_i = x;
        photo_CO₂_dependence!(photo_set, leaf);

        return FT(0)
    else
        # calculate flow rate and optimizer at the given p_i - Δp
        leaf.p_i = x - FT(1e-4);
        photo_CO₂_dependence!(photo_set, leaf);
        g_lc = leaf.An*FT(1e-6) / (p_a - leaf.p_i) * p_atm;
        g_sc = 1 / max( 1/g_lc - 1/g_m - 1/g_bc, FT(1e-3) );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        e_1  = g_lw * (p_sat - p_H₂O) / p_atm;
        k_1  = leaf_xylem_risk(hs, e_1);
        opt1 = kr_max * leaf.An + a_max * k_1;
    
        # calculate flow rate and optimizer at the given p_i
        leaf.p_i = x;
        photo_CO₂_dependence!(photo_set, leaf);
        g_lc = leaf.An*FT(1e-6) / (p_a - leaf.p_i) * p_atm;
        g_sc = 1 / max( 1/g_lc - 1/g_m - 1/g_bc, FT(1e-3) );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        e_2  = g_lw * (p_sat - p_H₂O) / p_atm;
        k_2  = leaf_xylem_risk(hs, e_2);
        opt2 = kr_max * leaf.An + a_max * k_2;
    
        # calculate the marginal gain - marginal penalty
        slope = (opt2 - opt1) / (e_2 - e_1);
    
        return slope
    end
end

function envir_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::OSMWang
            ) where {FT<:AbstractFloat}
    @unpack ec, g_bc, g_bw, g_m, p_sat = leaf;
    @unpack p_a, p_atm, p_H₂O = envir;

    # Wang model may run into numerical issues when E is too small
    # Thus, a control of maximal p_i = p_a applies just in case
    if (x >= p_a) || (isnan(x))
        x = p_a;
        leaf.p_i = x;
        photo_CO₂_dependence!(photo_set, leaf);

        return FT(0)
    else
        # calculate flow rate and optimizer at the given p_i - Δp
        leaf.p_i = x - FT(1e-4);
        photo_CO₂_dependence!(photo_set, leaf);
        g_lc = leaf.An*FT(1e-6) / (p_a - leaf.p_i) * p_atm;
        g_sc = 1 / max( 1/g_lc - 1/g_m - 1/g_bc, FT(1e-3) );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        e_1  = g_lw * (p_sat - p_H₂O) / p_atm;
        opt1 = leaf.An * (ec - e_1);
    
        # calculate flow rate and optimizer at the given p_i
        leaf.p_i = x;
        photo_CO₂_dependence!(photo_set, leaf);
        g_lc = leaf.An*FT(1e-6) / (p_a - leaf.p_i) * p_atm;
        g_sc = 1 / max( 1/g_lc - 1/g_m - 1/g_bc, FT(1e-3) );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        e_2  = g_lw * (p_sat - p_H₂O) / p_atm;
        opt2 = leaf.An * (ec - e_2);
    
        # calculate the marginal gain - marginal penalty
        slope = (opt2 - opt1) / (e_2 - e_1);
    
        return slope
    end
end

function envir_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::OSMWAP
            ) where {FT<:AbstractFloat}
    @unpack g_bc, g_bw, g_m, hs, p_sat = leaf;
    @unpack p_a, p_atm, p_H₂O = envir;

    # WAP model may run into numerical issues when E is too small
    # Thus, a control of maximal p_i = p_a applies just in case
    if (x >= p_a) || (isnan(x))
        x = p_a;
        leaf.p_i = x;
        photo_CO₂_dependence!(photo_set, leaf);

        return FT(0)
    else
        # calculate flow rate and optimizer at the given p_i - Δp
        leaf.p_i = x - FT(1e-4);
        photo_CO₂_dependence!(photo_set, leaf);
        g_lc = leaf.An*FT(1e-6) / (p_a - leaf.p_i) * p_atm;
        g_sc = 1 / max( 1/g_lc - 1/g_m - 1/g_bc, FT(1e-3) );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        e_1  = g_lw * (p_sat - p_H₂O) / p_atm;
        p_1  = xylem_p_from_flow(hs, e_1);
        opt1 = leaf.An - (sm.a * p_1^2 - sm.b * p_1);
    
        # calculate flow rate and optimizer at the given p_i
        leaf.p_i = x;
        photo_CO₂_dependence!(photo_set, leaf);
        g_lc = leaf.An*FT(1e-6) / (p_a - leaf.p_i) * p_atm;
        g_sc = 1 / max( 1/g_lc - 1/g_m - 1/g_bc, FT(1e-3) );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        e_2  = g_lw * (p_sat - p_H₂O) / p_atm;
        p_2  = xylem_p_from_flow(hs, e_1);
        opt2 = leaf.An - (sm.a * p_2^2 - sm.b * p_2);
    
        # calculate the marginal gain - marginal penalty
        slope = (opt2 - opt1) / (e_2 - e_1) * FT(1e-3);
    
        return slope
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
            sm::EmpiricalStomatalModel
            ) where {FT<:AbstractFloat}
    # update TD only if T changes, update p_sat as well
    if leaf.T != leaf.T_old
        photo_temperature_dependence!(photo_set, leaf, envir);
        leaf.T_old = leaf.T;
    end
    photo_radiation_dependence!(photo_set, leaf);

    # solve for optimal CO₂, A and g_sw updated here
    _sm = SecantMethod{FT}(leaf.Γ_star, leaf.Γ_star+2);
    _cs = CompactSolution();
    _st = SolutionTolerance{FT}(1e-2);
    @inline f(x) = envir_diff!(x, photo_set, leaf, envir, sm);
    _solut = find_zero(f, _sm, _cs, _st);

    # update leaf conductances
    leaf.g_lc = leaf.An*FT(1e-6) / (envir.p_a - leaf.p_i) * envir.p_atm;
    leaf.g_sc = 1 / max( 1/leaf.g_lc - 1/leaf.g_m - 1/leaf.g_bc, FT(1e-3) );
    leaf.g_sw = leaf.g_sc * FT(1.6);
    leaf.g_lw = 1 / (1/leaf.g_sw + 1/leaf.g_bw);
    
    # make sure g_sw in its range, and update flow rate
    leaf_gsw_control!(photo_set, leaf, envir);

    return nothing
end

function leaf_photo_from_envir!(
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            sm::OSMEller
            ) where {FT<:AbstractFloat}
    # update TD only if T changes, update p_sat as well
    if leaf.T != leaf.T_old
        photo_temperature_dependence!(photo_set, leaf, envir);
        leaf.T_old = leaf.T;
    end
    photo_radiation_dependence!(photo_set, leaf);

    # if An at p_a <= 0, close stomata; else, optimize it
    leaf.p_i = envir.p_a;
    photo_CO₂_dependence!(photo_set, leaf);
    if leaf.An > 0.05
        # define a starting point to ensure solution
        _pl = leaf.Γ_star;
        _dp = (envir.p_a - FT(0.01) - leaf.Γ_star) / 20;
        while true
            _pl     += _dp
            leaf.p_i = _pl;
            photo_CO₂_dependence!(photo_set, leaf);
            if leaf.An > 0
                break
            end
        end
        _pu = min(_pl+1, (envir.p_a+_pl)/2)

        # solve for optimal CO₂, A and g_sw updated here
        _sm = SecantMethod{FT}(_pl, _pu);
        _cs = CompactSolution();
        _st = SolutionTolerance{FT}(1e-2);
        @inline f(x) = envir_diff!(x, photo_set, leaf, envir, sm);
        _solut = find_zero(f, _sm, _cs, _st);

        leaf.g_lc = leaf.An*FT(1e-6) / (envir.p_a - leaf.p_i) * envir.p_atm;
        leaf.g_sc = 1 / max(1/leaf.g_lc - 1/leaf.g_m - 1/leaf.g_bc, FT(1e-3));
        leaf.g_sw = leaf.g_sc * FT(1.6);
        leaf.g_lw = 1 / (1/leaf.g_sw + 1/leaf.g_bw);
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
    _hs::LeafHydraulics = leaf.hs
    # if T changes, update TD and ec, then T_old
    if leaf.T != leaf.T_old
        photo_temperature_dependence!(photo_set, leaf, envir);
        leaf.ec = leaf_e_crit(_hs, leaf.ec);
        leaf.T_old = leaf.T;
    end
    photo_radiation_dependence!(photo_set, leaf);

    # if p_ups changes, update ec, then p_old
    if leaf.p_ups != leaf.p_old
        leaf.ec = leaf_e_crit(_hs, leaf.ec);
        leaf.p_old = leaf.p_ups;
    end

    # if An at p_a <= 0, close stomata; else, optimize it
    leaf.p_i = envir.p_a;
    photo_CO₂_dependence!(photo_set, leaf);
    if leaf.An > 0.05
        # reculate A_max and Kr_max whenever envir changes
        _g_crit     = leaf.ec / (leaf.p_sat - envir.p_H₂O) * envir.p_atm;
        _g_sw       = 1 / max(1/_g_crit - 1/leaf.g_bw, FT(1.6e-3));
        leaf.g_sw   = min(_g_crit, leaf.g_max);
        leaf.g_lw   = 1/ (1/leaf.g_bw + 1/leaf.g_sw);
        leaf.g_sc   = _g_sw / FT(1.6);
        leaf.g_lc   = 1 / (1/leaf.g_bc + 1/leaf.g_sc + 1/leaf.g_m);
        leaf_photo_from_glc!(photo_set, leaf, envir);
        leaf.a_max  = leaf.An;
        leaf.kr_max = _hs.k_history[end];

        # define a starting point to ensure solution
        _pl = leaf.Γ_star;
        _dp = (envir.p_a - FT(0.01) - leaf.Γ_star) / 20;
        while true
            _pl     += _dp
            leaf.p_i = _pl;
            photo_CO₂_dependence!(photo_set, leaf);
            if leaf.An > 0
                break
            end
        end
        _pu = min(_pl+1, (envir.p_a+_pl)/2)

        # solve for optimal CO₂, A and g_sw updated here
        _sm = SecantMethod{FT}(_pl, _pu);
        _cs = CompactSolution();
        _st = SolutionTolerance{FT}(1e-2);
        @inline f(x) = envir_diff!(x, photo_set, leaf, envir, sm);
        _solut = find_zero(f, _sm, _cs, _st);

        leaf.g_lc = leaf.An*FT(1e-6) / (envir.p_a - leaf.p_i) * envir.p_atm;
        leaf.g_sc = 1 / max(1/leaf.g_lc - 1/leaf.g_m - 1/leaf.g_bc, FT(1e-3));
        leaf.g_sw = leaf.g_sc * FT(1.6);
        leaf.g_lw = 1 / (1/leaf.g_sw + 1/leaf.g_bw);
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
            sm::OSMWang
            ) where {FT<:AbstractFloat}
    # if T changes, update TD and ec, then T_old
    if leaf.T != leaf.T_old
        photo_temperature_dependence!(photo_set, leaf, envir);
        leaf.ec = leaf_e_crit(leaf.hs, leaf.ec);
        leaf.T_old = leaf.T;
    end
    photo_radiation_dependence!(photo_set, leaf);

    # if p_ups changes, update ec, then p_old
    if leaf.p_ups != leaf.p_old
        leaf.ec = leaf_e_crit(leaf.hs, leaf.ec);
        leaf.p_old = leaf.p_ups;
    end

    # if An at p_a <= 0, close stomata; else, optimize it
    leaf.p_i = envir.p_a;
    photo_CO₂_dependence!(photo_set, leaf);
    if leaf.An > 0.05
        # define a starting point to ensure solution
        _pl = leaf.Γ_star;
        _dp = (envir.p_a - FT(0.01) - leaf.Γ_star) / 20;
        while true
            _pl     += _dp
            leaf.p_i = _pl;
            photo_CO₂_dependence!(photo_set, leaf);
            if leaf.An > 0
                break
            end
        end
        _pu = min(_pl+1, (envir.p_a+_pl)/2)

        # solve for optimal CO₂, A and g_sw updated here
        _sm = SecantMethod{FT}(_pl, _pu);
        _cs = CompactSolution();
        _st = SolutionTolerance{FT}(1e-2);
        @inline f(x) = envir_diff!(x, photo_set, leaf, envir, sm);
        _solut = find_zero(f, _sm, _cs, _st);

        leaf.g_lc = leaf.An*FT(1e-6) / (envir.p_a - leaf.p_i) * envir.p_atm;
        leaf.g_sc = 1 / max(1/leaf.g_lc - 1/leaf.g_m - 1/leaf.g_bc, FT(1e-3));
        leaf.g_sw = leaf.g_sc * FT(1.6);
        leaf.g_lw = 1 / (1/leaf.g_sw + 1/leaf.g_bw);
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
            sm::OSMWAP
            ) where {FT<:AbstractFloat}
    # update TD only if T changes, update p_sat as well
    if leaf.T != leaf.T_old
        photo_temperature_dependence!(photo_set, leaf, envir);
        leaf.T_old = leaf.T;
    end
    photo_radiation_dependence!(photo_set, leaf);

    # solve for optimal CO₂, A and g_sw updated here
    _sm = SecantMethod{FT}(leaf.Γ_star, leaf.Γ_star+2);
    _cs = CompactSolution();
    _st = SolutionTolerance{FT}(1e-2);
    @inline f(x) = envir_diff!(x, photo_set, leaf, envir, sm);
    _solut = find_zero(f, _sm, _cs, _st);

    # update leaf conductances
    leaf.g_lc = leaf.An*FT(1e-6) / (envir.p_a - leaf.p_i) * envir.p_atm;
    leaf.g_sc = 1 / max(1/leaf.g_lc - 1/leaf.g_m - 1/leaf.g_bc, FT(1e-3));
    leaf.g_sw = leaf.g_sc * FT(1.6);
    leaf.g_lw = 1 / (1/leaf.g_sw + 1/leaf.g_bw);
    
    # make sure g_sw in its range, and update flow rate
    leaf_gsw_control!(photo_set, leaf, envir);

    return nothing
end
