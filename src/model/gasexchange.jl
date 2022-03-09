###############################################################################
#
# Update CanopyLayer struct from stomatal conductance
#
###############################################################################
"""
    gas_exchange!(
                canopyi::CanopyLayer{FT},
                hs::LeafHydraulics{FT},
                psoil::FT,
                swc::FT,
                envir::AirLayer{FT},
                sm::EmpiricalStomatalModel{FT},
                bt::AbstractBetaFunction{FT}
    ) where {FT<:AbstractFloat}
    gas_exchange!(
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
                canopyi::CanopyLayer{FT},
                hs::LeafHydraulics{FT},
                envir::AirLayer{FT},
                sm::AbstractStomatalModel{FT}
    ) where {FT<:AbstractFloat}
    gas_exchange!(
                canopyi::CanopyLayer{FT},
                hs::LeafHydraulics{FT},
                envir::AirLayer{FT},
                sm::AbstractStomatalModel{FT},
                ind::Int
    ) where {FT<:AbstractFloat}
    gas_exchange!(
                canopyi::CanopyLayer{FT},
                hs::TreeSimple{FT},
                envir::AirLayer{FT},
                sm::OSMWang{FT}
    ) where {FT<:AbstractFloat}
    gas_exchange!(
                canopyi::CanopyLayer{FT},
                envir::AirLayer{FT},
                drive::GlcDrive,
                ind::Int,
                glc::FT
    ) where {FT<:AbstractFloat}
    gas_exchange!(
                canopyi::CanopyLayer{FT},
                envir::AirLayer{FT},
                drive::GlcDrive,
                ind::Int
    ) where {FT<:AbstractFloat}
    gas_exchange!(
                canopyi::CanopyLayer{FT},
                envir::AirLayer{FT},
                drive::GlcDrive
    ) where {FT<:AbstractFloat}
    gas_exchange!(
                canopyi::CanopyLayer{FT},
                envir::AirLayer{FT},
                drive::GswDrive,
                ind::Int,
                gsw::FT
    ) where {FT<:AbstractFloat}
    gas_exchange!(
                canopyi::CanopyLayer{FT},
                envir::AirLayer{FT},
                drive::GswDrive,
                ind::Int
    ) where {FT<:AbstractFloat}
    gas_exchange!(
                canopyi::CanopyLayer{FT},
                envir::AirLayer{FT},
                drive::GswDrive
    ) where {FT<:AbstractFloat}

Calculate steady state gas exchange rates, given
- `canopyi` [`CanopyLayer`](@ref) type struct
- `hs` Leaf hydraulic system or TreeSimple hydraulic organism
- `envir` [`AirLayer`] type struct
- `sm` [`EmpiricalStomatalModel`](@ref) or [`OptimizationStomatalModel`](@ref)
- `bt` [`AbstractBetaFunction`](@ref) type struct
- `ind` Nth leaf in canopyi
- `drive` [`GlcDrive`](@ref) or [`GswDrive`](@ref) drive mode
- `glc` Given leaf diffusive conductance to CO₂
- `gsw` Given stomatal conductance to H₂O

Note 1: When there is no drive mode in the parameter list, the function
    calculates the steady state stomatal conductance first, and then the gas
    exchange rates.

Note 2: When using GlcDrive mode, gas exchange rates are computed using the
    given `glc`. However, this option does not make the gsw control, so it is
    not guaranteed that gsw is within the physiological range. Thus, gsw
    control should be made outside this function. This option is supposed to be
    used in the optimal stomatl conductance models only, because optimal
    conductance can be outside the physiological stomatal conductance range.
    Thus, using this option for other purposes need to be cautious. In this
    case, it is recommended to use the GswDrive mode.

Note 3: When using GswDrive mode, gas exchange rates are computed using the
    given `gsw`. Moreover, `gsw` control so that gsw is within the
    physiological range.
"""
function gas_exchange!(
            canopyi::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            svc::AbstractSoilVC{FT},
            psoil::FT,
            swc::FT,
            envir::AirLayer{FT},
            sm::EmpiricalStomatalModel{FT},
            bt::AbstractBetaFunction{FT}
) where {FT<:AbstractFloat}
    # update the temperature dependent parameters
    update_leaf_TP!(canopyi, hs, envir);

    # calculate optimal solution for each leaf
    for ind in eachindex(canopyi.APAR)
        canopyi.ps.apar = canopyi.APAR[ind];
        gas_exchange!(canopyi, hs, svc, psoil, swc, envir, sm, bt, ind);
    end
end




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
            hs::LeafHydraulics{FT},
            envir::AirLayer{FT},
            sm::OptimizationStomatalModel{FT}
) where {FT<:AbstractFloat}
    # update the temperature dependent parameters and maximal a and kr
    update_leaf_TP!(canopyi, hs, envir);
    update_leaf_AK!(canopyi, hs, envir);

    # calculate optimal solution for each leaf
    for ind in eachindex(canopyi.APAR)
        canopyi.ps.apar = canopyi.APAR[ind];
        gas_exchange!(canopyi, hs, envir, sm, ind);
    end

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
            hs::TreeSimple{FT},
            envir::AirLayer{FT},
            sm::OSMWang{FT}
) where {FT<:AbstractFloat}
    # update the temperature dependent parameters and maximal a and kr
    update_leaf_TP!(canopyi, hs, envir);
    update_leaf_AK!(canopyi, hs.leaf, envir);

    # calculate optimal solution for each leaf
    for ind in eachindex(canopyi.APAR)
        canopyi.ps.apar = canopyi.APAR[ind];
        gas_exchange!(canopyi, hs.leaf, envir, sm, ind);
    end

    return nothing
end




function gas_exchange!(
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            drive::GlcDrive,
            ind::Int,
            glc::FT
) where {FT<:AbstractFloat}
    # update the conductances
    canopyi.g_lc[ind] = glc;
    gas_exchange!(canopyi, envir, drive, ind);

    return nothing
end




function gas_exchange!(
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            drive::GlcDrive,
            ind::Int
) where {FT<:AbstractFloat}
    # update the conductances
    canopyi.g_sc[ind] = 1 / ( 1 / canopyi.g_lc[ind] -
                              1 / canopyi.g_m[ind] -
                              1 / canopyi.g_bc[ind] );
    canopyi.g_sw[ind] = canopyi.g_sc[ind] * FT(1.6);
    canopyi.g_lw[ind] = 1 / ( 1 / canopyi.g_sw[ind] +
                              1 / canopyi.g_bw[ind] );

    # update the photosynthetic rates
    if canopyi.g_lc[ind] != canopyi.ps.g_CO₂
        leaf_photosynthesis!(canopyi.ps, envir, GCO₂Mode(), canopyi.g_lc[ind]);
    end
    canopyi.Ac[ind] = canopyi.ps.PSM.a_c;
    canopyi.Aj[ind] = canopyi.ps.PSM.a_j;
    canopyi.Ap[ind] = canopyi.ps.PSM.a_p;
    canopyi.Ag[ind] = canopyi.ps.PSM.a_gross;
    canopyi.An[ind] = canopyi.ps.PSM.a_net;
    canopyi.φs[ind] = canopyi.ps.PRC.ϕ_f;

    # update the pressures
    canopyi.p_i[ind] = canopyi.ps.p_CO₂_i;
    canopyi.p_s[ind] = canopyi.ps.p_CO₂_s;

    return nothing
end




function gas_exchange!(
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            drive::GlcDrive
) where {FT<:AbstractFloat}
    # update the conductances for each "leaf"
    for i in eachindex(canopyi.g_lc)
        canopyi.ps.apar = canopyi.APAR[i];
        gas_exchange!(canopyi, envir, drive, i, canopyi.g_lc[i]);
    end

    return nothing
end




function gas_exchange!(
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            drive::GswDrive,
            ind::Int,
            gsw::FT
) where {FT<:AbstractFloat}
    # update the conductances
    canopyi.g_sw[ind] = gsw;
    gas_exchange!(canopyi, envir, drive, ind);

    return nothing
end




function gas_exchange!(
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            drive::GswDrive,
            ind::Int
) where {FT<:AbstractFloat}
    @unpack g_max, g_min = canopyi;

    # update the conductances
    if canopyi.g_sw[ind] > g_max
        canopyi.g_sw[ind] = g_max;
    end

    if canopyi.g_sw[ind] < g_min
        canopyi.g_sw[ind] = g_min;
    end

    canopyi.g_lw[ind] = 1 / ( 1 / canopyi.g_sw[ind] +
                              1 / canopyi.g_bw[ind] );
    canopyi.g_sc[ind] = canopyi.g_sw[ind] / FT(1.6);
    canopyi.g_lc[ind] = 1 / ( FT(1.6) / canopyi.g_sw[ind] +
                              1 / canopyi.g_m[ind] +
                              1 / canopyi.g_bc[ind] );

    # update the photosynthetic rates
    leaf_photosynthesis!(canopyi.ps, envir, GCO₂Mode(), canopyi.g_lc[ind]);
    canopyi.Ac[ind] = canopyi.ps.PSM.a_c;
    canopyi.Aj[ind] = canopyi.ps.PSM.a_j;
    canopyi.Ap[ind] = canopyi.ps.PSM.a_p;
    canopyi.Ag[ind] = canopyi.ps.PSM.a_gross;
    canopyi.An[ind] = canopyi.ps.PSM.a_net;
    canopyi.φs[ind] = canopyi.ps.PRC.ϕ_f;

    # update the pressures
    canopyi.p_i[ind] = canopyi.ps.p_CO₂_i;
    canopyi.p_s[ind] = canopyi.ps.p_CO₂_s;

    return nothing
end




function gas_exchange!(
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            drive::GswDrive
) where {FT<:AbstractFloat}
    # update the conductances for each "leaf"
    for i in eachindex(canopyi.g_sw)
        canopyi.ps.apar = canopyi.APAR[i];
        gas_exchange!(canopyi, envir, drive, i);
    end

    return nothing
end
