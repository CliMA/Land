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
                svc::AbstractSoilVC{FT},
                psoil::FT,
                swc::FT,
                envir::AirLayer{FT},
                sm::OptimizationStomatalModel{FT},
                bt::AbstractBetaFunction{FT},
                mode::GlcDrive,
                ind::Int
    ) where {FT<:AbstractFloat}
    solution_diff!(x::FT,
                photo_set::AbstractPhotoModelParaSet{FT},
                canopyi::CanopyLayer{FT},
                hs::LeafHydraulics{FT},
                envir::AirLayer{FT},
                sm::AbstractStomatalModel{FT},
                mode::AbstractDrive,
                ind::Int
    ) where {FT<:AbstractFloat}

Calculate the difference to be minimized for a given
- `x` Assumed leaf diffusive conductance or stomatal conductance, depending on
    `mode`
- `photo_set`[`C3ParaSet`] or [`C4ParaSet`] type parameter set
- `canopyi`[`CanopyLayer`](@ref) type struct
- `hs` Leaf hydraulic system
- `psoil` Soil water potential `[MPa]`
- `swc` Soil water content
- `envir`[`AirLayer`] type struct
- `sm` [`EmpiricalStomatalModel`](@ref) or [`OptimizationStomatalModel`](@ref)
- `bt` [`AbstractBetaFunction`](@ref) type struct
- `mode` [`GlcDrive`](@ref) or [`GswDrive`](@ref) mode
- `ind` Nth leaf in the canopy layer

The former function works for all empirical stomatal models, and the latter
    works for all optimization based models.
"""
function solution_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            svc::AbstractSoilVC{FT},
            psoil::FT,
            swc::FT,
            envir::AirLayer{FT},
            sm::EmpiricalStomatalModel{FT},
            bt::AbstractBetaG{FT},
            mode::GlcDrive,
            ind::Int
) where {FT<:AbstractFloat}
    # unpack variables
    @unpack ps = canopyi;
    g_bc  = canopyi.g_bc[ind];
    g_m   = canopyi.g_m[ind];

    # update photosynthesis for ps
    leaf_photosynthesis!(photo_set, ps, envir, GCO₂Mode(), x);

    # calculate g_sw from stomatal model
    β    = β_factor(hs, svc, bt, hs.p_element[end], psoil, swc);
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
            svc::AbstractSoilVC{FT},
            psoil::FT,
            swc::FT,
            envir::AirLayer{FT},
            sm::EmpiricalStomatalModel{FT},
            bt::AbstractBetaV{FT},
            mode::GlcDrive,
            ind::Int
) where {FT<:AbstractFloat}
    # unpack variables
    ps   = canopyi.ps;
    g_bc = canopyi.g_bc[ind];
    g_m  = canopyi.g_m[ind];

    # update photosynthesis for ps
    leaf_photosynthesis!(photo_set, ps, envir, GCO₂Mode(), x);

    # make beta correction over the photosynthesis system
    β    = β_factor(hs, svc, bt, hs.p_element[end], psoil, swc);
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




function solution_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet,
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            envir::AirLayer{FT},
            sm::OSMEller,
            mode::GlcDrive,
            ind::Int
) where {FT<:AbstractFloat}
    # unpack variables
    @unpack g_max, g_min, p_sat, ps = canopyi;
    @unpack p_atm, p_H₂O = envir;
    g_bc = canopyi.g_bc[ind];
    g_bw = canopyi.g_bw[ind];
    g_m  = canopyi.g_m[ind];

    # calculate the limit of g_lc to ensure x in the physiological range
    _glh = 1 / (1/g_bc + FT(1.6)/g_max + 1/g_m);
    _gll = 1 / (1/g_bc + FT(1.6)/g_min + 1/g_m);

    if x > _glh
        x = _glh + FT(1e-4);
        canopyi.g_sw[ind] = 2 * g_max;

        return FT(0)
    elseif x< _gll
        x = _gll - FT(1e-4);
        canopyi.g_sw[ind] = FT(0);

        return FT(0)
    else
        # update photosynthesis from x-FT(1e-3)
        g_lc = x - FT(1e-3);
        leaf_photosynthesis!(photo_set, ps, envir, GCO₂Mode(), g_lc);

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        a_1  = ps.An;
        e_1  = g_lw * (p_sat - p_H₂O) / p_atm;
        k_1  = xylem_risk(hs, e_1);

        # update photosynthesis from x
        g_lc = x;
        leaf_photosynthesis!(photo_set, ps, envir, GCO₂Mode(), g_lc);

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        a_2  = ps.An;
        e_2  = g_lw * (p_sat - p_H₂O) / p_atm;
        k_2  = xylem_risk(hs, e_2);

        ∂A∂E = (a_2 - a_1) / (e_2 - e_1);
        ∂Θ∂E = (k_1 - k_2) / (e_2 - e_1) * a_2 / k_2;
        diff = ∂A∂E - ∂Θ∂E;

        #= used for debugging
        @show x;
        @show ∂A∂E;
        @show ∂Θ∂E;
        println();
        #sleep(0.1);
        # =#

        return diff
    end
end




function solution_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            envir::AirLayer{FT},
            sm::OSMSperry{FT},
            mode::GlcDrive,
            ind::Int
) where {FT<:AbstractFloat}
    # unpack variables
    @unpack g_max, g_min, kr_max, p_sat, ps = canopyi;
    @unpack p_atm, p_H₂O = envir;
    a_max = canopyi.a_max[ind];
    g_bc  = canopyi.g_bc[ind];
    g_bw  = canopyi.g_bw[ind];
    g_m   = canopyi.g_m[ind];

    # calculate the limit of g_lc to ensure x in the physiological range
    _glh = 1 / (1/g_bc + FT(1.6)/g_max + 1/g_m);
    _gll = 1 / (1/g_bc + FT(1.6)/g_min + 1/g_m);

    if x > _glh
        x = _glh + FT(1e-4);
        canopyi.g_sw[ind] = 2 * g_max;

        return FT(0)
    elseif x< _gll
        x = _gll - FT(1e-4);
        canopyi.g_sw[ind] = FT(0);

        return FT(0)
    else
        # update photosynthesis from x-FT(1e-3)
        g_lc = x - FT(1e-3);
        leaf_photosynthesis!(photo_set, ps, envir, GCO₂Mode(), g_lc);

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        a_1  = ps.An;
        e_1  = g_lw * (p_sat - p_H₂O) / p_atm;
        k_1  = xylem_risk(hs, e_1);

        # update photosynthesis from x
        g_lc = x;
        leaf_photosynthesis!(photo_set, ps, envir, GCO₂Mode(), g_lc);

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        a_2  = ps.An;
        e_2  = g_lw * (p_sat - p_H₂O) / p_atm;
        k_2  = xylem_risk(hs, e_2);

        ∂A∂E = (a_2 - a_1) / (e_2 - e_1);
        ∂Θ∂E = (k_1 - k_2) / (e_2 - e_1) * a_max / kr_max;
        diff = ∂A∂E - ∂Θ∂E;

        #= used for debugging
        @show x;
        @show ∂A∂E;
        @show ∂Θ∂E;
        println();
        #sleep(0.1);
        # =#

        return diff
    end
end




function solution_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            envir::AirLayer{FT},
            sm::OSMWang{FT},
            mode::GlcDrive,
            ind::Int
) where {FT<:AbstractFloat}
    # unpack variables
    @unpack ec, g_max, g_min, p_sat, ps = canopyi;
    @unpack p_atm, p_H₂O = envir;
    g_bc  = canopyi.g_bc[ind];
    g_bw  = canopyi.g_bw[ind];
    g_m   = canopyi.g_m[ind];

    # calculate the limit of g_lc to ensure x in the physiological range
    _glh = 1 / (1/g_bc + FT(1.6)/g_max + 1/g_m);
    _gll = 1 / (1/g_bc + FT(1.6)/g_min + 1/g_m);

    if x > _glh
        x = _glh + FT(1e-4);
        canopyi.g_sw[ind] = 2 * g_max;

        return FT(0)
    elseif x< _gll
        x = _gll - FT(1e-4);
        canopyi.g_sw[ind] = FT(0);

        return FT(0)
    else
        # update photosynthesis from x-FT(1e-3)
        g_lc = x - FT(1e-3);
        leaf_photosynthesis!(photo_set, ps, envir, GCO₂Mode(), g_lc);

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        a_1  = ps.An;
        e_1  = g_lw * (p_sat - p_H₂O) / p_atm;

        # update photosynthesis from x
        g_lc = x;
        leaf_photosynthesis!(photo_set, ps, envir, GCO₂Mode(), g_lc);

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        a_2  = ps.An;
        e_2  = g_lw * (p_sat - p_H₂O) / p_atm;

        ∂A∂E = (a_2 - a_1) / (e_2 - e_1);
        ∂Θ∂E = a_2 / max(ec - e_2, FT(1e-7));
        diff = ∂A∂E - ∂Θ∂E;

        #= used for debugging
        @show x;
        @show ∂A∂E;
        @show ∂Θ∂E;
        println();
        #sleep(0.1);
        # =#

        return diff
    end
end




function solution_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            envir::AirLayer{FT},
            sm::OSMWAP{FT},
            mode::GlcDrive,
            ind::Int
) where {FT<:AbstractFloat}
    # unpack variables
    @unpack g_max, g_min, p_sat, ps = canopyi;
    @unpack p_atm, p_H₂O = envir;
    g_bc = canopyi.g_bc[ind];
    g_bw = canopyi.g_bw[ind];
    g_m  = canopyi.g_m[ind];

    # calculate the limit of g_lc to ensure x in the physiological range
    _glh = 1 / (1/g_bc + FT(1.6)/g_max + 1/g_m);
    _gll = 1 / (1/g_bc + FT(1.6)/g_min + 1/g_m);

    if x > _glh
        x = _glh + FT(1e-4);
        canopyi.g_sw[ind] = 2 * g_max;

        return FT(0)
    elseif x< _gll
        x = _gll - FT(1e-4);
        canopyi.g_sw[ind] = FT(0);

        return FT(0)
    else
        # update photosynthesis from x-FT(1e-3)
        g_lc = x - FT(1e-3);
        leaf_photosynthesis!(photo_set, ps, envir, GCO₂Mode(), g_lc);

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        a_1  = ps.An;
        e_1  = g_lw * (p_sat - p_H₂O) / p_atm;
        p_1  = end_pressure(hs, e_1);

        # update photosynthesis from x
        g_lc = x;
        leaf_photosynthesis!(photo_set, ps, envir, GCO₂Mode(), g_lc);

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        a_2  = ps.An;
        e_2  = g_lw * (p_sat - p_H₂O) / p_atm;
        p_2  = end_pressure(hs, e_2);

        ∂A∂E = (a_2 - a_1) / (e_2 - e_1);
        ∂Θ∂E = (-2 * sm.a * p_2 + sm.b) / (e_2 - e_1) * (p_1 - p_2);
        diff = ∂A∂E - ∂Θ∂E;

        #= used for debugging
        @show x;
        @show ∂A∂E;
        @show ∂Θ∂E;
        println();
        #sleep(0.1);
        # =#

        return diff
    end
end




function solution_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            envir::AirLayer{FT},
            sm::OSMWAPMod{FT},
            mode::GlcDrive,
            ind::Int
) where {FT<:AbstractFloat}
    # unpack variables
    @unpack g_max, g_min, p_sat, ps = canopyi;
    @unpack p_atm, p_H₂O = envir;
    g_bc = canopyi.g_bc[ind];
    g_bw = canopyi.g_bw[ind];
    g_m  = canopyi.g_m[ind];

    # calculate the limit of g_lc to ensure x in the physiological range
    _glh = 1 / (1/g_bc + FT(1.6)/g_max + 1/g_m);
    _gll = 1 / (1/g_bc + FT(1.6)/g_min + 1/g_m);

    if x > _glh
        x = _glh + FT(1e-4);
        canopyi.g_sw[ind] = 2 * g_max;

        return FT(0)
    elseif x< _gll
        x = _gll - FT(1e-4);
        canopyi.g_sw[ind] = FT(0);

        return FT(0)
    else
        # update photosynthesis from x-FT(1e-3)
        g_lc = x - FT(1e-3);
        leaf_photosynthesis!(photo_set, ps, envir, GCO₂Mode(), g_lc);

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        a_1  = ps.An;
        e_1  = g_lw * (p_sat - p_H₂O) / p_atm;
        p_1  = end_pressure(hs, e_1);

        # update photosynthesis from x
        g_lc = x;
        leaf_photosynthesis!(photo_set, ps, envir, GCO₂Mode(), g_lc);

        g_sc = 1 / ( 1/g_lc - 1/g_m - 1/g_bc );
        g_sw = g_sc * FT(1.6);
        g_lw = 1 / (1/g_sw + 1/g_bw);
        a_2  = ps.An;
        e_2  = g_lw * (p_sat - p_H₂O) / p_atm;
        p_2  = end_pressure(hs, e_2);

        ∂A∂E = (a_2 - a_1) / (e_2 - e_1);
        ∂Θ∂E = (-2 * sm.a * p_2 * a_2) / (e_2 - e_1) * (p_1 - p_2);
        diff = ∂A∂E - ∂Θ∂E;

        #= used for debugging
        @show x;
        @show ∂A∂E;
        @show ∂Θ∂E;
        println();
        #sleep(0.1);
        # =#

        return diff
    end
end




function solution_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            envir::AirLayer{FT},
            sm::OptimizationStomatalModel{FT},
            mode::GswDrive,
            ind::Int
) where {FT<:AbstractFloat}
    # gsw has already been guaranteed to be in the range
    # unpack variables
    @unpack ec, g_max, g_min, p_sat, ps = canopyi;
    @unpack p_atm, p_H₂O = envir;
    g_bc = canopyi.g_bc[ind];
    g_bw = canopyi.g_bw[ind];
    g_m  = canopyi.g_m[ind];
    g_lc = 1 / (1/g_bc + FT(1.6)/x + 1/g_m);

    return solution_diff!(g_lc, photo_set, canopyi, hs, envir, sm, GlcDrive(),
                          ind)
end
