# These functions are used to find steady state solution, move them to SPAC.jl
#=
function gas_exchange!(
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            svc::AbstractSoilVC{FT},
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
        @unpack P_AIR, p_H₂O = envir;
        _g_bc = canopyi.g_bc[ind];
        _g_bw = canopyi.g_bw[ind];
        _g_m  = canopyi.g_m[ind];

        # solve for optimal g_lc, A and g_sw updated here
        _gh = 1 / (1/_g_bc + FT(1.6)/g_max + 1/_g_m);
        _gl = 1 / (1/_g_bc + FT(1.6)/g_min + 1/_g_m);
        _sm = NewtonBisectionMethod{FT}(x_min=_gl, x_max=_gh);
        _st = SolutionTolerance{FT}(1e-4, 50);
        @inline f(x) = solution_diff!(x, canopyi, hs, svc, psoil, swc, envir, sm, bt, GlcDrive(), ind);
        _solut = find_zero(f, _sm, _st);

        # update leaf conductances and rates
        gas_exchange!(canopyi, envir, GlcDrive(), ind, _solut);
    else
        canopyi.g_sw[ind] = FT(0);
    end

    # make sure g_sw in its range
    gsw_control!(canopyi, envir, ind);

    return nothing
end



function gas_exchange!(
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            sm::EmpiricalStomatalModel{FT},
            ind::Int
) where {FT<:AbstractFloat}
    if canopyi.APAR[ind] > 1
        # unpack required variables
        @unpack ec, g_max, g_min, p_sat = canopyi;
        @unpack P_AIR, p_H₂O = envir;
        _g_bc = canopyi.g_bc[ind];
        _g_bw = canopyi.g_bw[ind];
        _g_m  = canopyi.g_m[ind];

        # solve for optimal g_lc, A and g_sw updated here
        _gh = 1 / (1/_g_bc + FT(1.6)/g_max + 1/_g_m);
        _gl = 1 / (1/_g_bc + FT(1.6)/g_min + 1/_g_m);
        _sm = NewtonBisectionMethod{FT}(x_min=_gl, x_max=_gh);
        _st = SolutionTolerance{FT}(1e-4, 50);
        @inline f(x) = solution_diff!(x, canopyi, envir, sm, ind);
        _solut = find_zero(f, _sm, _st);

        # update leaf conductances and rates
        gas_exchange!(canopyi, envir, GlcDrive(), ind, _solut);
    else
        canopyi.g_sw[ind] = FT(0);
    end

    # make sure g_sw in its range
    gsw_control!(canopyi, envir, ind);

    return nothing
end



function gas_exchange!(
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            envir::AirLayer{FT},
            sm::OptimizationStomatalModel{FT},
            ind::Int
) where {FT<:AbstractFloat}
    # if there is light
    if canopyi.APAR[ind] > 1
        # unpack required variables
        @unpack ec, g_max, g_min, p_sat = canopyi;
        @unpack P_AIR, p_H₂O = envir;
        _g_bc   = canopyi.g_bc[ind];
        _g_bw   = canopyi.g_bw[ind];
        _g_m    = canopyi.g_m[ind];

        # calculate the physiological maximal g_sw
        _g_crit = ec / (p_sat - p_H₂O) * P_AIR;
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
            _gh = 1 / (1/_g_bc + FT(1.6)/_g_max + 1/_g_m);
            _gl = 1 / (1/_g_bc + FT(1.6)/g_min  + 1/_g_m);
            _sm = NewtonBisectionMethod{FT}(x_min=_gl, x_max=_gh);
            _st = SolutionTolerance{FT}(1e-4, 50);
            @inline fd(x) = solution_diff!(x, canopyi, hs, envir, sm, GlcDrive(), ind);
            _sl = find_zero(fd, _sm, _st);

            #= used for debugging
            @show _g_sw;
            @show _g_max;
            @show leaf.p_ups;
            @show _solut;
            println();
            =#

            # update leaf conductances and rates
            gas_exchange!(canopyi, envir, GlcDrive(), ind, _sl);
        end

    # if there is no light, use nighttime mode
    else
        @unpack ec, g_max, g_min, p_sat = canopyi;
        @unpack P_AIR, p_H₂O = envir;
        _sm = NewtonBisectionMethod{FT}(x_min=g_min, x_max=g_max);
        _st = SolutionTolerance{FT}(1e-4, 50);
        @inline fn(x) = nocturnal_diff!(x, canopyi, envir, sm);
        _sl = find_zero(fn, _sm, _st);

        # update leaf conductances and rates
        gas_exchange!(canopyi, envir, GswDrive(), ind, _sl);
    end

    # make sure g_sw in its range
    gsw_control!(canopyi, envir, ind);

    return nothing
end




function gas_exchange!(
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
        @unpack P_AIR, p_H₂O = envir;
        _g_bc   = canopyi.g_bc[ind];
        _g_bw   = canopyi.g_bw[ind];
        _g_m    = canopyi.g_m[ind];

        # calculate the physiological maximal g_sw
        # add an 90% offset in _g_crit for Eller model to avoid dK/dE = 0
        _g_crit = FT(0.9) * ec / (p_sat - p_H₂O) * P_AIR;
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
            _gh = 1 / (1/_g_bc + FT(1.6)/_g_max + 1/_g_m);
            _gl = 1 / (1/_g_bc + FT(1.6)/g_min  + 1/_g_m);
            _sm = NewtonBisectionMethod{FT}(x_min=_gl, x_max=_gh);
            _st = SolutionTolerance{FT}(1e-4, 50);
            @inline f(x) = solution_diff!(x, canopyi, hs, envir, sm, GlcDrive(), ind);
            _solut = find_zero(f, _sm, _st);

            #= used for debugging
            @show _g_sw;
            @show _g_max;
            @show leaf.p_ups;
            @show _solut;
            println();
            =#

            # update leaf conductances and rates
            gas_exchange!(canopyi, envir, GlcDrive(), ind, _solut);
        end

    # if there is no light
    else
        canopyi.g_sw[ind] = FT(0);
    end

    # make sure g_sw in its range
    gsw_control!(canopyi, envir, ind);

    return nothing
end




function gas_exchange!(
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            envir::AirLayer{FT},
            sm::OSMSperry{FT},
            ind::Int
) where {FT<:AbstractFloat}
    # if there is light
    if canopyi.APAR[ind] > 1
        # unpack required variables
        @unpack ec, g_max, g_min, p_sat = canopyi;
        @unpack P_AIR, p_H₂O = envir;
        _g_bc   = canopyi.g_bc[ind];
        _g_bw   = canopyi.g_bw[ind];
        _g_m    = canopyi.g_m[ind];

        # calculate the physiological maximal g_sw
        # add an 70% offset in _g_crit for Eller model to avoid dK/dE = 0
        _g_crit = FT(0.7) * ec / (p_sat - p_H₂O) * P_AIR;
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
            _gh = 1 / (1/_g_bc + FT(1.6)/_g_max + 1/_g_m);
            _gl = 1 / (1/_g_bc + FT(1.6)/g_min  + 1/_g_m);
            _sm = NewtonBisectionMethod{FT}(x_min=_gl, x_max=_gh);
            _st = SolutionTolerance{FT}(1e-4, 50);
            @inline f(x) = solution_diff!(x, canopyi, hs, envir, sm, GlcDrive(), ind);
            _solut = find_zero(f, _sm, _st);

            #= used for debugging
            @show _g_sw;
            @show _g_max;
            @show leaf.p_ups;
            @show _solut;
            println();
            =#

            # update leaf conductances and rates
            gas_exchange!(canopyi, envir, GlcDrive(), ind, _solut);
        end

    # if there is no light
    else
        canopyi.g_sw[ind] = FT(0);
    end

    # make sure g_sw in its range
    gsw_control!(canopyi, envir, ind);

    return nothing
end




function gas_exchange!(
            canopyi::CanopyLayer{FT},
            spac::MonoElementSPAC{FT},
            envir::AirLayer{FT},
            sm::OSMWang{FT}
) where {FT<:AbstractFloat}
    # update the temperature dependent parameters and maximal a and kr
    update_leaf_TP!(canopyi, spac, envir);
    update_leaf_AK!(canopyi, spac.LEAF.HS, envir);

    # calculate optimal solution for each leaf
    for ind in eachindex(canopyi.APAR)
        canopyi.ps.apar = canopyi.APAR[ind];
        gas_exchange!(canopyi, spac.LEAF.HS, envir, sm, ind);
    end

    return nothing
end
=#
