###############################################################################
#
# Initialize SPACMono structure
#
###############################################################################
"""
    initialize_spac_canopy!(node::SPACMono{FT}) where {FT<:AbstractFloat}

Initialize the RT parameters for a given
- `node` [`SPACMono`](@ref) type struct
"""
function initialize_spac_canopy!(node::SPACMono{FT}) where {FT<:AbstractFloat}
    # 0.1 create variables required
    @unpack angles, envirs, in_rad, leaves_rt, n_canopy, plant_ps, photo_set, rt_con, soil_opt, wl_set = node;
    canopy_rt = node.canopy_rt;
    can_opt = node.can_opt;
    can_rad = node.can_rad;
    plant_hs = node.plant_hs;
    fraction_sl::Array{FT,1} = repeat(canopy_rt.lidf, outer=[canopy_rt.nAzi]) / length(canopy_rt.lazitab);
    n_sl = length(canopy_rt.lidf) * length(canopy_rt.lazitab);

    # fluspect the canopy layers
    for ican in 1:n_canopy
        fluspect!(leaves_rt[ican], wl_set);
    end

    # Four Different steps to compute Short-Wave RT
    canopy_geometry!(canopy_rt, angles, can_opt, rt_con)
    canopy_matrices!(leaves_rt, can_opt);
    short_wave!(canopy_rt, can_opt, can_rad, in_rad, soil_opt, rt_con);
    canopy_fluxes!(canopy_rt, can_opt, can_rad, in_rad, soil_opt, leaves_rt, wl_set, rt_con);

    # Compute Long Wave (Last term is LW incoming in W m^-2)
    thermal_fluxes!(leaves_rt, can_opt, can_rad, canopy_rt, soil_opt, [FT(400.0)], wl_set);

    # update the canopy leaf area partition information
    for i_can in 1:n_canopy
        rt_layer = n_canopy + 1 - i_can;
        iPS      = plant_ps[i_can];

        # calculate the fraction of sunlit and shaded leaves
        f_view = (can_opt.Ps[rt_layer] + can_opt.Ps[rt_layer+1]) / 2;

        for iLF in 1:n_sl
            iPS.APAR[iLF] = can_rad.absPAR_sunCab[(rt_layer-1)*n_sl + iLF] * FT(1e6);
            iPS.LAIx[iLF] = f_view * fraction_sl[iLF];
        end
        iPS.APAR[end] = can_rad.absPAR_shadeCab[rt_layer] * FT(1e6);
        iPS.LAIx[end] = 1 - f_view;

        update_leaf_TP!(photo_set, iPS, plant_hs.leaves[i_can], envirs[i_can]);

        envirs[i_can].t_air = T₂₅(FT);
        envirs[i_can].p_sat = saturation_vapor_pressure(T₂₅(FT));
        envirs[i_can].p_H₂O = envirs[i_can].p_sat / 2;
    end

    return nothing
end








###############################################################################
#
# As the land model is modularized, there could be multiple places that conatin
#     exact the same information. Thus, these information needs to be synced.
#     This type of information includes leaf area, hydraulic traits, and
#     photosynthetic traits.
#
# The following function is meant for leaf area
#
###############################################################################
"""
    update_LAI!(node::SPACMono{FT}, lai::FT) where {FT<:AbstractFloat}

Update SPAC leaf area index
"""
function update_LAI!(node::SPACMono{FT}, lai::FT) where {FT<:AbstractFloat}
    # 1. Update the general information of the SPAC
    node.la = node.ga * lai;

    # 2. Update the LAI for canopy radiation model
    node.canopy_rt.LAI = lai;

    # 3. Update the LA and LAI for stomatal and photosynthesis models
    for _iPS in node.plant_ps
        _iPS.LA   = node.la / node.n_canopy;
        _iPS.LAI  = lai / node.n_canopy;
        _iPS.tLAI = lai;
    end

    # 4. Update the LA for hydraulic model
    #    TODO abstract this to avoid memory allocation
    for _iHS in node.plant_hs.leaves
        _iHS.area = node.la / node.n_canopy;
    end

    return nothing
end




"""
    update_Kmax!(node::SPACMono{FT}, kmax::FT) where {FT<:AbstractFloat}

Update the maximal hydraulic conductance for SPAC

TODO need to abstractize this to PlantHydraulics.jl
"""
function update_Kmax!(node::SPACMono{FT}, kmax::FT) where {FT<:AbstractFloat}
    # Root:Stem:Leaf conductance ratio is set as 1:2:2.
    #    For trees, roots: 2X, trunk: 8X, branch: 8X, leaves: 4X
    #    For palms, roots: 2X, trunk: 4X, leaves: 4X
    #    For grasses, roots: 2X, leaves: 2X

    # 1. update the hydraulic conductance in Roots
    for root in node.plant_hs.roots
        root.k_max      = 2 * kmax / node.plant_hs.n_root;
        root.k_element .= root.k_max * root.N;
    end

    # 2. If is a tree
    if typeof(node.plant_hs) <: TreeLikeOrganism
        node.plant_hs.trunk.k_max      = 8 * kmax;
        node.plant_hs.trunk.k_element .= node.plant_hs.trunk.k_max * node.plant_hs.trunk.N;
        for stem in node.plant_hs.branch
            stem.k_max      = 8 * kmax / node.plant_hs.n_canopy;
            stem.k_element .= stem.k_max * stem.N;
        end
        for leaf in node.plant_hs.leaves
            leaf.k_sla      = 4 * kmax / node.la;
            leaf.k_element .= leaf.k_sla * leaf.N;
        end
    end

    # 3. If is a palm
    if typeof(node.plant_hs) <: PalmLikeOrganism
        node.plant_hs.trunk.k_max      = 4 * kmax;
        node.plant_hs.trunk.k_element .= node.plant_hs.trunk.k_max * node.plant_hs.trunk.N;
        for leaf in node.plant_hs.leaves
            leaf.k_sla      = 4 * kmax / node.la;
            leaf.k_element .= leaf.k_sla * leaf.N;
        end
    end

    # 4. If is a grass
    if typeof(node.plant_hs) <: GrassLikeOrganism
        for leaf in node.plant_hs.leaves
            leaf.k_sla      = 2 * kmax / node.la;
            leaf.k_element .= leaf.k_sla * leaf.N;
        end
    end

    return nothing
end




"""
    update_VJRWW!(node::SPACMono{FT}, vcmax::FT) where {FT<:AbstractFloat}

Update Vcmax25(WW), Jmax25(WW), and Rd25(WW) from a given Vcmax25
"""
function update_VJRWW!(node::SPACMono{FT}, vcmax::FT; expo::FT = FT(NaN)) where {FT<:AbstractFloat}
    # TODO change the ratio accordingly to photo_set
    # TODO add another ratio V2J in photo_set
    # Update Vcmax25, Jmax25 (1.67 Vcmax), and Rd25 (0.015 Vcmax)
    for _iPS in node.plant_ps
        _iPS.ps.Vcmax     = vcmax;
        _iPS.ps.Vcmax25   = vcmax;
        _iPS.ps.Vcmax25WW = vcmax;
        _iPS.ps.Vpmax     = vcmax;
        _iPS.ps.Vpmax25   = vcmax;
        _iPS.ps.Vpmax25WW = vcmax;
        _iPS.ps.Jmax      = vcmax * 1.67;
        _iPS.ps.Jmax25    = vcmax * 1.67;
        _iPS.ps.Jmax25WW  = vcmax * 1.67;
        _iPS.ps.Rd        = vcmax * 0.015;
        _iPS.ps.Rd25      = vcmax * 0.015;
        _iPS.ps.Rd25WW    = vcmax * 0.015;
    end

    if isnan(expo)
        return nothing
    end;

    # if expo is not NaN, add a exponential trend
    for _i_can in eachindex(node.plant_ps)
        _x_rt = exp( -expo * node.canopy_rt.LAI * (1 - _i_can / node.n_canopy) );
        node.plant_ps[_i_can].ps.Vcmax25WW *= _x_rt;
        node.plant_ps[_i_can].ps.Vpmax25WW *= _x_rt;
        node.plant_ps[_i_can].ps.Jmax25WW  *= _x_rt;
        node.plant_ps[_i_can].ps.Rd25WW    *= _x_rt;
        node.plant_ps[_i_can].ps.Vcmax25   *= _x_rt;
        node.plant_ps[_i_can].ps.Vpmax25   *= _x_rt;
        node.plant_ps[_i_can].ps.Jmax25    *= _x_rt;
        node.plant_ps[_i_can].ps.Rd25      *= _x_rt;
        node.plant_ps[_i_can].ps.Vcmax     *= _x_rt;
        node.plant_ps[_i_can].ps.Vpmax     *= _x_rt;
        node.plant_ps[_i_can].ps.Jmax      *= _x_rt;
        node.plant_ps[_i_can].ps.Rd        *= _x_rt;
    end;

    return nothing
end




"""
    update_VJR!(node::SPACMono{FT}, ratio::FT) where {FT<:AbstractFloat}

Update effective Vcmax25, Jmax25, and Rd25 from a given tune factor
"""
function update_VJR!(node::SPACMono{FT}, ratio::FT) where {FT<:AbstractFloat}
    # TODO change the ratio accordingly to photo_set
    # TODO add another ratio V2J in photo_set
    # Update Vcmax25, Jmax25 (1.67 Vcmax), and Rd25 (0.015 Vcmax)
    for _iPS in node.plant_ps
        _iPS.ps.Vcmax25 = ratio * _iPS.ps.Vcmax25WW;
        _iPS.ps.Jmax25  = ratio * _iPS.ps.Jmax25WW;
        _iPS.ps.Rd25    = ratio * _iPS.ps.Rd25WW;
        _iPS.ps.Vcmax   = ratio * _iPS.ps.Vcmax25WW;
        _iPS.ps.Jmax    = ratio * _iPS.ps.Jmax25WW;
        _iPS.ps.Rd      = ratio * _iPS.ps.Rd25WW;
    end

    return nothing
end




"""
    update_Cab!(node::SPACMono{FT}, cab::FT) where {FT<:AbstractFloat}

Update leaf chlorophyll content, and then rerun `fluspect!`, given
- `cab` chlorophyll content
- `car` carotenoid content
"""
function update_Cab!(node::SPACMono{FT}, cab::FT; cab_2_car::FT = FT(1/7)) where {FT<:AbstractFloat}
    for _leaf in node.leaves_rt
        _leaf.Cab = cab;
        _leaf.Car = cab * cab_2_car;
        fluspect!(_leaf, node.wl_set);
    end

    return nothing
end




"""
    update_Weibull!(node::SPACMono{FT}, b::FT, c::FT) where {FT<:AbstractFloat}

Update Weibull B and C for all components
"""
function update_Weibull!(node::SPACMono{FT}, b::FT) where {FT<:AbstractFloat}
    # 1. update B and C for roots
    for _root in node.plant_hs.roots
        _root.vc.b = b;
    end

    # 2. update B and C for trunk (if exists)
    if !(typeof(node.plant_hs) <: GrassLikeOrganism)
        node.plant_hs.trunk.vc.b = b;
    end

    # 3. update B and C for branches (if exist)
    if typeof(node.plant_hs) <: TreeLikeOrganism
        for _stem in node.plant_hs.branch
            _stem.vc.b = b;
        end
    end

    # 4. update B and C for leaves
    for _leaf in node.plant_hs.leaves
        _leaf.vc.b = b;
    end

    return nothing
end




function update_Weibull!(node::SPACMono{FT}, b::FT, c::FT) where {FT<:AbstractFloat}
    # 1. update B and C for roots
    for _root in node.plant_hs.roots
        _root.vc.b = b;
        _root.vc.c = c;
    end

    # 2. update B and C for trunk (if exists)
    if !(typeof(node.plant_hs) <: GrassLikeOrganism)
        node.plant_hs.trunk.vc.b = b;
        node.plant_hs.trunk.vc.c = c;
    end

    # 3. update B and C for branches (if exist)
    if typeof(node.plant_hs) <: TreeLikeOrganism
        for _stem in node.plant_hs.branch
            _stem.vc.b = b;
            _stem.vc.c = c;
        end
    end

    # 4. update B and C for leaves
    for _leaf in node.plant_hs.leaves
        _leaf.vc.b = b;
        _leaf.vc.c = c;
    end

    return nothing
end


"""
    sync_par!(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Sync canopy layers PAR and sunlit fractions
"""
function sync_par!(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    # calculate leaf level flux per canopy layer
    _nSL = spac.canopy_rt.nAzi * spac.canopy_rt.nIncl;
    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can];
        _iRT = spac.n_canopy + 1 - _i_can;

        # calculate the fraction of sunlit and shaded leaves
        _f_view = (spac.can_opt.Ps[_iRT] + spac.can_opt.Ps[_iRT+1]) / 2;
        for iLF in 1:_nSL
            _iPS.APAR[iLF] = spac.can_rad.absPAR_sunCab[(_iRT-1)*_nSL+iLF] * FT(1e6);
            _iPS.LAIx[iLF] = _f_view * spac.f_SL[iLF];
        end;
        _iPS.APAR[end] = spac.can_rad.absPAR_shadeCab[_iRT] * FT(1e6);
        _iPS.LAIx[end] = 1 - _f_view;
    end;

    return nothing
end


"""
    sync_fqy!(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Sync canopy layer fluorescence quantum yield
"""
function sync_fqy!(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    # update fluorescence quantum yield for all modes
    for _i_can in 1:spac.n_canopy
        _iRT = spac.n_canopy + 1 - _i_can;
        _iPS = spac.plant_ps[_i_can];
        spac.can_rad.ϕ_sun[:,:,_iRT] .= reshape(view(_iPS.φs,1:spac.canopy_rt.nIncl*spac.canopy_rt.nAzi), spac.canopy_rt.nIncl, spac.canopy_rt.nAzi);
        spac.can_rad.ϕ_shade[_iRT] = (_iPS).φs[end];
    end

    return nothing
end


"""
    update_par!(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Compute PAR and sync it to SPAC
"""
function update_par!(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    canopy_geometry!(spac.canopy_rt, spac.angles, spac.can_opt, spac.rt_con);
    canopy_matrices!(spac.leaves_rt, spac.can_opt);
    short_wave!(spac.canopy_rt, spac.can_opt, spac.can_rad, spac.in_rad, spac.soil_opt, spac.rt_con);
    canopy_fluxes!(spac.canopy_rt, spac.can_opt, spac.can_rad, spac.in_rad, spac.soil_opt, spac.leaves_rt, spac.wl_set, spac.rt_con);
    sync_par!(spac);

    return nothing
end


"""
    update_sif!(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Compute SIF and sync it to SPAC
"""
function update_sif!(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    sync_fqy!(spac);
    SIF_fluxes!(spac.leaves_rt, spac.can_opt, spac.can_rad, spac.canopy_rt, spac.soil_opt, spac.wl_set, spac.rt_con, spac.rt_dim);

    return nothing
end


"""
    prescribe_air!(spac::SPACMono{FT}, co2::FT, p_atm::FT, t_air::FT, vpd::FT, wind::FT) where {FT<:AbstractFloat}

Prescribe environmental conditions, given
- `spac` SPAC
- `co2` CO2 in `[ppm]`
- `p_atm` Atmospheric pressure in `[Pa]`
- `t_air` Air temperature in `[K]`
- `vpd` Atmospheric vapor pressure deficit
- `wind` Wind speed in `[m s⁻¹]`

"""
function prescribe_air!(spac::SPACMono{FT}, co2::FT, p_atm::FT, t_air::FT, vpd::FT, wind::FT) where {FT<:AbstractFloat}
    for _i_can in 1:spac.n_canopy
        _iEN = spac.envirs[_i_can];

        # update environmental conditions
        _iEN.t_air = t_air;
        _iEN.p_atm = p_atm;
        _iEN.p_a   = _iEN.p_atm * co2 * 1e-6;
        _iEN.p_O₂  = _iEN.p_atm * 0.209;
        _iEN.p_sat = saturation_vapor_pressure(_iEN.t_air);
        _iEN.vpd   = vpd;
        _iEN.p_H₂O = _iEN.p_sat - _iEN.vpd;
        _iEN.RH    = _iEN.p_H₂O / _iEN.p_sat;
        _iEN.wind  = wind;
    end;

    return nothing
end


"""
    prescribe_t_leaf!(spac::SPACMono{FT}, t_leaf::FT) where {FT<:AbstractFloat}

Prescribe leaf temperature, given
- `spac` SPAC
- `t_leaf` Leaf temperature in `[K]`

"""
function prescribe_t_leaf!(spac::SPACMono{FT}, t_leaf::FT) where {FT<:AbstractFloat}
    for _i_can in 1:spac.n_canopy
        _iEN = spac.envirs[_i_can];
        _iHS = spac.plant_hs.leaves[_i_can];
        _iPS = spac.plant_ps[_i_can];

        # prescribe leaf temperature
        _iPS.T = t_leaf;
        update_leaf_TP!(spac.photo_set, _iPS, _iHS, _iEN);
        temperature_effects!(_iHS, _iPS.T);
    end;

    return nothing
end


"""
    prescribe_swc!(spac::SPACMono{FT}, args...) where {FT<:AbstractFloat}

Prescribe soil moisture, given
- `spac` SPAC
- `args` A list of soil moisture (from 0 to Θs)

"""
function prescribe_swc!(spac::SPACMono{FT}, args...) where {FT<:AbstractFloat}
    for _i in eachindex(spac.plant_hs.roots)
        _svc = spac.plant_hs.roots[_i].sh;
        spac.swc[_i] = max(_svc.Θr + eps(FT), args[_i]);
        spac.plant_hs.roots[_i].p_ups = soil_p_25_swc(_svc, spac.swc[_i]);
    end

    return nothing
end
