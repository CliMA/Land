###############################################################################
#
# Diff function to minimize by ConstrainedRootSolvers
# Useful for Ballberry, Leuning, and Medlyn models
#
###############################################################################
"""
    solution_diff!(x::FT,
                photo_set::AbstractPhotoModelParaSet{FT},
                canopyi::CanopyLayer{FT},
                hs::LeafHydraulics{FT},
                psoil::FT,
                swc::FT,
                envir::AirLayer{FT},
                sm::OptimizationStomatalModel{FT},
                bt::AbstractBetaFunction{FT},
                ind::Int
    ) where {FT<:AbstractFloat}
    solution_diff!(x::FT,
                photo_set::AbstractPhotoModelParaSet{FT},
                canopyi::CanopyLayer{FT},
                hs::LeafHydraulics{FT},
                envir::AirLayer{FT},
                sm::AbstractStomatalModel{FT},
                ind::Int
    ) where {FT<:AbstractFloat}

Calculate the difference to be minimized for a given
- `x` Assumed leaf diffusive conductance
- `photo_set`[`C3ParaSet`] or [`C4ParaSet`] type parameter set
- `canopyi`[`CanopyLayer`](@ref) type struct
- `hs` Leaf hydraulic system
- `psoil` Soil water potential `[MPa]`
- `swc` Soil water content
- `envir`[`AirLayer`] type struct
- `sm` [`EmpiricalStomatalModel`](@ref) or [`OptimizationStomatalModel`](@ref)
- `bt` [`AbstractBetaFunction`](@ref) type struct
- `ind` Nth leaf in the canopy layer

The former function works for Ball-Berry, Leuning, and Medlyn models, the
    latter works for Gentine and all the optimization based models.
"""
function solution_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            psoil::FT,
            swc::FT,
            envir::AirLayer{FT},
            sm::EmpiricalStomatalModel{FT},
            bt::AbstractBetaG{FT},
            ind::Int
) where {FT<:AbstractFloat}
    # unpack variables
    @unpack ps = canopyi;
    g_bc  = canopyi.g_bc[ind];
    g_m   = canopyi.g_m[ind];

    # update photosynthesis for ps
    leaf_photosynthesis!(photo_set, ps, envir, GCO₂Mode(), x);

    # calculate g_sw from stomatal model
    β    = β_factor(bt, hs.p_element[end], psoil, swc);
    g_md = stomatal_conductance(sm, ps, envir, β);
    g_md = min(canopyi.g_max, g_md);

    # calculate model predicted g_lc
    g_lm = 1 / (FT(1.6)/g_md + 1/g_bc + 1/g_m);

    return g_lm - x
end




function solution_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            psoil::FT,
            swc::FT,
            envir::AirLayer{FT},
            sm::EmpiricalStomatalModel{FT},
            bt::AbstractBetaV{FT},
            ind::Int
) where {FT<:AbstractFloat}
    # unpack variables
    ps   = canopyi.ps;
    g_bc = canopyi.g_bc[ind];
    g_m  = canopyi.g_m[ind];

    # update photosynthesis for ps
    leaf_photosynthesis!(photo_set, ps, envir, GCO₂Mode(), x);

    # make beta correction over the photosynthesis system
    β    = β_factor(bt, hs.p_element[end], psoil, swc);
    _rat = ps.Vcmax25WW * β / ps.Vcmax25;
    if _rat != 1
        ps.Jmax25  *= _rat;
        ps.Jmax    *= _rat;
        ps.Rd25    *= _rat;
        ps.Rd      *= _rat;
        ps.Vcmax25 *= _rat;
        ps.Vcmax   *= _rat;
        ps.Vpmax25 *= _rat;
        ps.Vpmax   *= _rat;
    end

    # calculate g_sw from stomatal model
    g_md = stomatal_conductance(sm, ps, envir, FT(1));
    g_md = min(canopyi.g_max, g_md);

    # calculate model predicted g_lc
    g_lm = 1 / (FT(1.6)/g_md + 1/g_bc + 1/g_m);

    return g_lm - x
end








###############################################################################
#
# Solve for stomatal conductance from environmental conditions
# Useful for BallBerry, Gentine, Leuning, and Medlyn models
#
###############################################################################
"""
    gas_exchange!(
                photo_set::AbstractPhotoModelParaSet{FT},
                canopyi::CanopyLayer{FT},
                hs::LeafHydraulics{FT},
                psoil::FT,
                swc::FT,
                envir::AirLayer{FT},
                sm::EmpiricalStomatalModel{FT},
                bt::AbstractBetaFunction{FT}
    ) where {FT<:AbstractFloat}
    gas_exchange!(
                photo_set::AbstractPhotoModelParaSet{FT},
                canopyi::CanopyLayer{FT},
                hs::LeafHydraulics{FT},
                psoil::FT,
                swc::FT,
                envir::AirLayer{FT},
                sm::EmpiricalStomatalModel{FT},
                bt::AbstractBetaFunction{FT},
                ind::Int
    ) where {FT<:AbstractFloat}
    gas_exchange!(
                photo_set::AbstractPhotoModelParaSet{FT},
                canopyi::CanopyLayer{FT},
                hs::LeafHydraulics{FT},
                envir::AirLayer{FT},
                sm::AbstractStomatalModel{FT}
    ) where {FT<:AbstractFloat}
    gas_exchange!(
                photo_set::AbstractPhotoModelParaSet{FT},
                canopyi::CanopyLayer{FT},
                hs::LeafHydraulics{FT},
                envir::AirLayer{FT},
                sm::AbstractStomatalModel{FT},
                ind::Int
    ) where {FT<:AbstractFloat}

Calculate steady state gsw and photosynthesis from empirical approach, given
- `photo_set` [`C3ParaSet`] or [`C4ParaSet`] type parameter set
- `canopyi` [`CanopyLayer`](@ref) type struct
- `hs` Leaf hydraulic system
- `envir` [`AirLayer`] type struct
- `sm` [`EmpiricalStomatalModel`](@ref) or [`OptimizationStomatalModel`](@ref)
- `bt` [`AbstractBetaFunction`](@ref) type struct
- `ind` Nth leaf in canopyi

The former function works for Ball-Berry, Leuning, and Medlyn models, the
    latter works for Gentine and all the optimization based models.
"""
function gas_exchange!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            psoil::FT,
            swc::FT,
            envir::AirLayer{FT},
            sm::EmpiricalStomatalModel{FT},
            bt::AbstractBetaFunction{FT}
) where {FT<:AbstractFloat}
    # update the temperature dependent parameters
    update_leaf_TP!(photo_set, canopyi, hs, envir);

    # calculate optimal solution for each leaf
    for ind in eachindex(canopyi.APAR)
        canopyi.ps.APAR = canopyi.APAR[ind];
        leaf_ETR!(photo_set, canopyi.ps);
        gas_exchange!(photo_set, canopyi, hs, psoil, swc, envir, sm, bt, ind);
    end
end




function gas_exchange!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            psoil::FT,
            swc::FT,
            envir::AirLayer{FT},
            sm::EmpiricalStomatalModel{FT},
            bt::AbstractBetaFunction{FT},
            ind::Int
) where {FT<:AbstractFloat}
    if canopyi.APAR[ind] > 1
        # unpack required variables
        @unpack ec, g_max, g_min, p_sat = canopyi;
        @unpack p_atm, p_H₂O = envir;
        _g_bc = canopyi.g_bc[ind];
        _g_bw = canopyi.g_bw[ind];
        _g_m  = canopyi.g_m[ind];

        # solve for optimal g_lc, A and g_sw updated here
        _gh    = 1 / (1/_g_bc + FT(1.6)/g_max + 1/_g_m);
        _gl    = 1 / (1/_g_bc + FT(1.6)/g_min + 1/_g_m);
        _sm    = NewtonBisectionMethod{FT}(_gl, _gh, (_gl+_gh)/2);
        _st    = SolutionTolerance{FT}(1e-4, 50);
        @inline f(x) = solution_diff!(x, photo_set, canopyi, hs, psoil, swc,
                                      envir, sm, bt, ind);
        _solut = find_zero(f, _sm, _st);

        # update leaf conductances and rates
        update_leaf_from_glc!(photo_set, canopyi, envir, ind, _solut);
    else
        canopyi.g_sw[ind] = FT(0);
    end

    # make sure g_sw in its range
    leaf_gsw_control!(photo_set, canopyi, envir, ind);

    return nothing
end
