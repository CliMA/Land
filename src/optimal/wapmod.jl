###############################################################################
#
# Diff function to minimize by ConstrainedRootSolvers
# Description in general model section in the empirical folder
#
###############################################################################
function solution_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            envir::AirLayer{FT},
            sm::OSMWAPMod{FT},
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
        println("");
        #sleep(0.1);
        # =#

        return diff
    end
end
