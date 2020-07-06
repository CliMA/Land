###############################################################################
#
# Solve for stomatal conductance from environmental conditions for all leaves
#
###############################################################################
"""
    leaf_photo_from_envir!(
            photo_set::AbstractPhotoModelParaSet,
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            sm::AbstractStomatalModel)
    leaf_photo_from_envir!(
            photo_set::AbstractPhotoModelParaSet,
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            sm::AbstractStomatalModel,
            ind::Int)

Calculate steady state gsw and photosynthesis from empirical approach, given
- `photo_set` [`C3ParaSet`] or [`C4ParaSet`] type parameter set
- `leaves` [`Leaves`](@ref) type struct
- `envir` [`AirLayer`] type struct
- `sm` [`EmpiricalStomatalModel`](@ref) or [`OptimizationStomatalModel`](@ref)
- `ind` Nth leaf in leaves
"""
function leaf_photo_from_envir!(
            photo_set::AbstractPhotoModelParaSet,
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            sm::OptimizationStomatalModel
            ) where {FT<:AbstractFloat}
    # update the temperature dependent parameters and maximal a and kr
    update_leaf_TP!(photo_set, leaves, envir);
    update_leaf_AK!(photo_set, leaves, envir);

    # calculate optimal solution for each leaf
    for ind in eachindex(leaves.APAR)
        leaf_photo_from_envir!(photo_set, leaves, envir, sm, ind);
    end
end

function leaf_photo_from_envir!(
            photo_set::AbstractPhotoModelParaSet,
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            sm::EmpiricalStomatalModel
            ) where {FT<:AbstractFloat}
    # update the temperature dependent parameters
    update_leaf_TP!(photo_set, leaves, envir);

    # calculate optimal solution for each leaf
    for ind in eachindex(leaves.APAR)
        leaf_photo_from_envir!(photo_set, leaves, envir, sm, ind);
    end
end

function leaf_photo_from_envir!(
            photo_set::AbstractPhotoModelParaSet,
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            sm::OSMEller,
            ind::Int
            ) where {FT<:AbstractFloat}
    # if there is light
    if leaves.APAR[ind] > 1
        # unpack required variables
        @unpack ec, g_max, g_min, p_sat = leaves;
        @unpack p_atm, p_H₂O = envir;
        _g_bc   = leaves.g_bc[ind];
        _g_bw   = leaves.g_bw[ind];
        _g_m    = leaves.g_m[ind];

        # calculate the physiological maximal g_sw
        # add an 90% offset in _g_crit for Eller model to avoid dK/dE = 0
        _g_crit = FT(0.9) * ec / (p_sat - p_H₂O) * p_atm;
        _g_max  = 1 / max(1/_g_crit - 1/_g_bw, FT(1e-3));
        _g_max  = min(_g_max, g_max);

        # if _g_sw is lower than g_min
        if _g_max <= g_min
            leaves.g_sw[ind] = FT(0);

        # if _g_sw is higher than g_min
        elseif leaves.a_max[ind] < 0.05
                leaves.g_sw[ind] = FT(0);
        else
            # solve for optimal g_lc, A and g_sw updated here
            _gh    = 1 / (1/_g_bc + FT(1.6)/_g_max + 1/_g_m);
            _gl    = 1 / (1/_g_bc + FT(1.6)/g_min  + 1/_g_m);
            _sm    = NewtonBisectionMethod{FT}(_gl, _gh, (_gl+_gh)/2);
            _st    = SolutionTolerance{FT}(1e-4);
            @inline f(x) = envir_diff!(x, photo_set, leaves, envir, sm, ind);
            _solut = find_zero(f, _sm, _st);

            #= used for debugging
            @show _g_sw;
            @show _g_max;
            @show leaf.p_ups;
            @show _solut;
            println("");
            =#

            # update leaf conductances and rates
            update_leaf_from_glc!(photo_set, leaves, envir, ind, _solut);
        end

    # if there is no light
    else
        leaves.g_sw[ind] = FT(0);
    end

    # make sure g_sw in its range
    leaf_gsw_control!(photo_set, leaves, envir, ind);

    return nothing
end

function leaf_photo_from_envir!(
            photo_set::AbstractPhotoModelParaSet,
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            sm::OSMSperry,
            ind::Int
            ) where {FT<:AbstractFloat}
    # if there is light
    if leaves.APAR[ind] > 1
        # unpack required variables
        @unpack ec, g_max, g_min, p_sat = leaves;
        @unpack p_atm, p_H₂O = envir;
        _g_bc   = leaves.g_bc[ind];
        _g_bw   = leaves.g_bw[ind];
        _g_m    = leaves.g_m[ind];

        # calculate the physiological maximal g_sw
        # add an 70% offset in _g_crit for Eller model to avoid dK/dE = 0
        _g_crit = FT(0.7) * ec / (p_sat - p_H₂O) * p_atm;
        _g_max  = 1 / max(1/_g_crit - 1/_g_bw, FT(1e-3));
        _g_max  = min(_g_max, g_max);

        # if _g_sw is lower than g_min
        if _g_max <= g_min
            leaves.g_sw[ind] = FT(0);

        # if _g_sw is higher than g_min
        elseif leaves.a_max[ind] < 0.05
                leaves.g_sw[ind] = FT(0);
        else
            # solve for optimal g_lc, A and g_sw updated here
            _gh    = 1 / (1/_g_bc + FT(1.6)/_g_max + 1/_g_m);
            _gl    = 1 / (1/_g_bc + FT(1.6)/g_min  + 1/_g_m);
            _sm    = NewtonBisectionMethod{FT}(_gl, _gh, (_gl+_gh)/2);
            _st    = SolutionTolerance{FT}(1e-4);
            @inline f(x) = envir_diff!(x, photo_set, leaves, envir, sm, ind);
            _solut = find_zero(f, _sm, _st);

            #= used for debugging
            @show _g_sw;
            @show _g_max;
            @show leaf.p_ups;
            @show _solut;
            println("");
            =#

            # update leaf conductances and rates
            update_leaf_from_glc!(photo_set, leaves, envir, ind, _solut);
        end

    # if there is no light
    else
        leaves.g_sw[ind] = FT(0);
    end

    # make sure g_sw in its range
    leaf_gsw_control!(photo_set, leaves, envir, ind);

    return nothing
end

function leaf_photo_from_envir!(
            photo_set::AbstractPhotoModelParaSet,
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            sm::OptimizationStomatalModel,
            ind::Int
            ) where {FT<:AbstractFloat}
    # if there is light
    if leaves.APAR[ind] > 1
        # unpack required variables
        @unpack ec, g_max, g_min, p_sat = leaves;
        @unpack p_atm, p_H₂O = envir;
        _g_bc   = leaves.g_bc[ind];
        _g_bw   = leaves.g_bw[ind];
        _g_m    = leaves.g_m[ind];

        # calculate the physiological maximal g_sw
        _g_crit = ec / (p_sat - p_H₂O) * p_atm;
        _g_max  = 1 / max(1/_g_crit - 1/_g_bw, FT(1e-3));
        _g_max  = min(_g_max, g_max);

        # if _g_sw is lower than g_min
        if _g_max <= g_min
            leaves.g_sw[ind] = FT(0);

        # if _g_sw is higher than g_min
        elseif leaves.a_max[ind] < 0.05
                leaves.g_sw[ind] = FT(0);
        else
            # solve for optimal g_lc, A and g_sw updated here
            _gh    = 1 / (1/_g_bc + FT(1.6)/_g_max + 1/_g_m);
            _gl    = 1 / (1/_g_bc + FT(1.6)/g_min  + 1/_g_m);
            _sm    = NewtonBisectionMethod{FT}(_gl, _gh, (_gl+_gh)/2);
            _st    = SolutionTolerance{FT}(1e-4);
            @inline f(x) = envir_diff!(x, photo_set, leaves, envir, sm, ind);
            _solut = find_zero(f, _sm, _st);

            #= used for debugging
            @show _g_sw;
            @show _g_max;
            @show leaf.p_ups;
            @show _solut;
            println("");
            =#

            # update leaf conductances and rates
            update_leaf_from_glc!(photo_set, leaves, envir, ind, _solut);
        end

    # if there is no light
    else
        leaves.g_sw[ind] = FT(0);
    end

    # make sure g_sw in its range
    leaf_gsw_control!(photo_set, leaves, envir, ind);

    return nothing
end

function leaf_photo_from_envir!(
            photo_set::AbstractPhotoModelParaSet,
            leaves::Leaves{FT},
            envir::AirLayer{FT},
            sm::EmpiricalStomatalModel,
            ind::Int
            ) where {FT<:AbstractFloat}
    if leaves.APAR[ind] > 1
        # unpack required variables
        @unpack ec, g_max, g_min, p_sat = leaves;
        @unpack p_atm, p_H₂O = envir;
        _g_bc   = leaves.g_bc[ind];
        _g_bw   = leaves.g_bw[ind];
        _g_m    = leaves.g_m[ind];

        # solve for optimal g_lc, A and g_sw updated here
        _gh    = 1 / (1/_g_bc + FT(1.6)/g_max + 1/_g_m);
        _gl    = 1 / (1/_g_bc + FT(1.6)/g_min + 1/_g_m);
        _sm    = NewtonBisectionMethod{FT}(_gl, _gh, (_gl+_gh)/2);
        _st    = SolutionTolerance{FT}(1e-4);
        @inline f(x) = envir_diff!(x, photo_set, leaves, envir, sm, ind);
        _solut = find_zero(f, _sm, _st);

        # update leaf conductances and rates
        update_leaf_from_glc!(photo_set, leaves, envir, ind, _solut);
    else
        leaves.g_sw[ind] = FT(0);
    end

    # make sure g_sw in its range
    leaf_gsw_control!(photo_set, leaves, envir, ind);

    return nothing
end
