###############################################################################
#
# Diff function to minimize by ConstrainedRootSolvers
# Description in BallBerry model section
#
###############################################################################
function solution_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet,
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            envir::AirLayer{FT},
            sm::OSMEller,
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
        println("");
        #sleep(0.1);
        # =#

        return diff
    end
end








###############################################################################
#
# Solve for stomatal conductance from environmental conditions
# Description in general model section in the empirical folder
#
###############################################################################
function gas_exchange!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            envir::AirLayer{FT},
            sm::OSMEller{FT},
            ind::Int
) where {FT<:AbstractFloat}
    # if there is light
    if canopyi.APAR[ind] > 1
        # unpack required variables
        @unpack ec, g_max, g_min, p_sat = canopyi;
        @unpack p_atm, p_H₂O = envir;
        _g_bc   = canopyi.g_bc[ind];
        _g_bw   = canopyi.g_bw[ind];
        _g_m    = canopyi.g_m[ind];

        # calculate the physiological maximal g_sw
        # add an 90% offset in _g_crit for Eller model to avoid dK/dE = 0
        _g_crit = FT(0.9) * ec / (p_sat - p_H₂O) * p_atm;
        _g_max  = 1 / max(1/_g_crit - 1/_g_bw, FT(1e-3));
        _g_max  = min(_g_max, g_max);

        # if _g_sw is lower than g_min
        if _g_max <= g_min
            canopyi.g_sw[ind] = FT(0);

        # if _g_sw is higher than g_min
        elseif canopyi.a_max[ind] < 0.05
                canopyi.g_sw[ind] = FT(0);
        else
            # solve for optimal g_lc, A and g_sw updated here
            _gh    = 1 / (1/_g_bc + FT(1.6)/_g_max + 1/_g_m);
            _gl    = 1 / (1/_g_bc + FT(1.6)/g_min  + 1/_g_m);
            _sm    = NewtonBisectionMethod{FT}(_gl, _gh, (_gl+_gh)/2);
            _st    = SolutionTolerance{FT}(1e-4, 50);
            @inline f(x) = solution_diff!(x, photo_set, canopyi, hs, envir, sm,
                                          ind);
            _solut = find_zero(f, _sm, _st);

            #= used for debugging
            @show _g_sw;
            @show _g_max;
            @show leaf.p_ups;
            @show _solut;
            println("");
            =#

            # update leaf conductances and rates
            update_leaf_from_glc!(photo_set, canopyi, envir, ind, _solut);
        end

    # if there is no light
    else
        canopyi.g_sw[ind] = FT(0);
    end

    # make sure g_sw in its range
    leaf_gsw_control!(photo_set, canopyi, envir, ind);

    return nothing
end
