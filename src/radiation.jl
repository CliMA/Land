#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-09: migrate the function from CanopyLayers
#     2022-Jun-09: rename function to canopy_radiation!
#     2022-Jun-10: add documentation
#
#######################################################################################################################################################################################################
"""
This function updates canopy radiation profiles. The supported methods include

$(METHODLIST)

"""
function canopy_radiation! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-09: migrate the function from CanopyLayers
#     2022-Jun-09: clean the function
#     2022-Jun-10: rename PAR/APAR to APAR/PPAR to be more accurate
#     2022-Jun-10: add PAR calculation (before absorption)
#     2022-Jun-10: add documentation
#     2022-Jun-10: compute shortwave net radiation
#
#######################################################################################################################################################################################################
"""

    canopy_radiation!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::HyperspectralRadiation{FT}, soil::Soil{FT}) where {FT<:AbstractFloat}

Updates canopy radiation profiles for shortwave radiation, given
- `can` `HyperspectralMLCanopy` type struct
- `leaves` Vector of `Leaf`
- `rad` Incoming solar radiation
- `soil` Bottom soil boundary layer
"""
canopy_radiation!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::HyperspectralRadiation{FT}, soil::Soil{FT}) where {FT<:AbstractFloat} = (
    @unpack APAR_CAR, N_LAYER, OPTICS, P_INCL, RADIATION, WLSET = can;
    _ilai = can.lai * can.ci / N_LAYER;
    _tlai = can.lai / N_LAYER;

    # 1. update upward and downward direct and diffuse radiation profiles
    RADIATION.e_direct[:,1] .= rad.e_direct;
    RADIATION.e_diffuse_down[:,1] .= rad.e_diffuse;

    for _i in 1:N_LAYER
        _e_d_i = view(RADIATION.e_diffuse_down,:,_i  );     # downward diffuse radiation at upper boundary
        _e_d_j = view(RADIATION.e_diffuse_down,:,_i+1);     # downward diffuse radiation at lower boundary
        _e_s_i = view(RADIATION.e_direct      ,:,_i  );     # direct radiation at upper boundary
        _e_s_j = view(RADIATION.e_direct      ,:,_i+1);     # direct radiation at lower boundary
        _e_u_i = view(RADIATION.e_diffuse_up  ,:,_i  );     # upward diffuse radiation at upper boundary

        _r_dd_i = view(OPTICS.ρ_dd,:,_i);   # reflectance of the upper boundary (i)
        _r_sd_i = view(OPTICS.ρ_sd,:,_i);   # reflectance of the upper boundary (i)
        _t_dd_i = view(OPTICS.τ_dd,:,_i);   # transmittance of the layer (i)
        _t_sd_i = view(OPTICS.τ_sd,:,_i);   # transmittance of the layer (i)
        _t_ss__ = OPTICS._τ_ss;             # transmittance for directional->directional

        _e_s_j .= _t_ss__ .* _e_s_i;
        _e_d_j .= _t_sd_i .* _e_s_i .+ _t_dd_i .* _e_d_i;
        _e_u_i .= _r_sd_i .* _e_s_i .+ _r_dd_i .* _e_d_i;
    end;

    _end = lastindex(OPTICS.ρ_sd, 2);
    RADIATION.e_diffuse_up[:,end] = view(OPTICS.ρ_sd,:,_end) .* view(RADIATION.e_direct,:,_end) .+ view(OPTICS.ρ_dd,:,_end) .* view(RADIATION.e_diffuse_down,:,_end);

    # 2. update the sunlit and shaded sum radiation and total absorbed radiation per layer and for soil
    for _i in 1:N_LAYER
        _a_s_i = view(RADIATION.e_net_direct  ,:,_i  );     # net absorbed direct radiation
        _a_d_i = view(RADIATION.e_net_diffuse ,:,_i  );     # net absorbed diffuse radiation
        _e_d_i = view(RADIATION.e_diffuse_down,:,_i  );     # downward diffuse radiation at upper boundary
        _e_s_i = view(RADIATION.e_direct      ,:,_i  );     # direct radiation at upper boundary
        _e_u_j = view(RADIATION.e_diffuse_up  ,:,_i+1);     # upward diffuse radiation at lower boundary
        _p_s_i = view(RADIATION.e_sum_direct  ,:,_i  );     # sum direct radiation
        _p_d_i = view(RADIATION.e_sum_diffuse ,:,_i  );     # sum diffuse radiation

        _r_dd_i = view(OPTICS.ρ_dd,:,_i);   # reflectance of the upper boundary (i)
        _r_sd_i = view(OPTICS.ρ_sd,:,_i);   # reflectance of the upper boundary (i)
        _t_dd_i = view(OPTICS.τ_dd,:,_i);   # transmittance of the layer (i)
        _t_sd_i = view(OPTICS.τ_sd,:,_i);   # transmittance of the layer (i)
        _t_ss__ = OPTICS._τ_ss;             # transmittance for directional->directional

        _p_s_i .= _e_s_i;
        _p_d_i .= _e_d_i .+ _e_u_j;

        _a_s_i .= _p_s_i .* (1 .- _t_ss__ .- _t_sd_i .- _r_sd_i);
        _a_d_i .= _p_d_i .* (1 .- _t_dd_i .- _r_dd_i);
    end;

    soil.e_net_direct .= view(RADIATION.e_direct,:,_end) .* (1 .- soil.ρ_sw);
    soil.e_net_diffuse .= view(RADIATION.e_diffuse_down,:,_end) .* (1 .- soil.ρ_sw);

    # 3. compute the spectra at the observer direction
    for _i in 1:N_LAYER
        _e_d_i = view(RADIATION.e_diffuse_down,:,_i);   # downward diffuse radiation at upper boundary
        _e_u_i = view(RADIATION.e_diffuse_up  ,:,_i);   # upward diffuse radiation at upper boundary

        _dob_i = view(OPTICS.σ_dob,:,_i);   # scattering coefficient backward for diffuse->observer
        _dof_i = view(OPTICS.σ_dob,:,_i);   # scattering coefficient forward for diffuse->observer
        _so__i = view(OPTICS.σ_so ,:,_i);   # bidirectional from solar to observer

        RADIATION.e_v[:,_i] .= (OPTICS.po[_i] .* _dob_i .* _e_d_i .+ OPTICS.po[_i] .* _dof_i .* _e_u_i .+ OPTICS.pso[_i] .* _so__i .* rad.e_direct) * _ilai;
    end;
    RADIATION.e_v[:,end] .= OPTICS.po[end] .* view(RADIATION.e_diffuse_up,:,_end);

    for _i in eachindex(RADIATION.e_o)
        RADIATION.e_o[_i] = sum(view(RADIATION.e_o,_i,:)) / FT(pi);
    end;

    RADIATION.albedo .= RADIATION.e_o * FT(pi) ./ (rad.e_direct .+ rad.e_diffuse);

    # 4. compute net absorption for leaves and soil
    for _i in 1:N_LAYER
        _Σ_shaded = view(RADIATION.e_net_diffuse,:,_i)' * WLSET.ΔΛ / 1000 / _tlai;
        _Σ_sunlit = view(RADIATION.e_net_direct ,:,_i)' * WLSET.ΔΛ / 1000 / _tlai;
        RADIATION.r_net_sw_shaded[_i] = _Σ_shaded;
        RADIATION.r_net_sw_sunlit[_i] = _Σ_sunlit / OPTICS.p_sunlit[_i] + _Σ_shaded;
        RADIATION.r_net_sw[_i] = _Σ_shaded + _Σ_sunlit;
    end;

    soil.r_net_sw = (soil.e_net_direct' * WLSET.ΔΛ + soil.e_net_diffuse' * WLSET.ΔΛ) / 1000;

    # 5. compute top-of-canopy and leaf level PAR, APAR, and PPAR
    RADIATION._par_shaded .= photon.(WLSET.Λ_PAR, view(rad.e_diffuse,WLSET.IΛ_PAR)) .* 1000;
    RADIATION._par_sunlit .= photon.(WLSET.Λ_PAR, view(rad.e_direct ,WLSET.IΛ_PAR)) .* 1000;
    RADIATION.par_in_diffuse = RADIATION._par_shaded' * WLSET.ΔΛ_PAR;
    RADIATION.par_in_direct = RADIATION._par_sunlit' * WLSET.ΔΛ_PAR;
    RADIATION.par_in = RADIATION.par_in_diffuse + RADIATION.par_in_direct;

    mul!(OPTICS._tmp_vec_azi, OPTICS._abs_fs', P_INCL);
    _normi = 1 / mean(OPTICS._tmp_vec_azi);

    for _i in 1:N_LAYER
        _α_apar = APAR_CAR ? view(leaves[_i].BIO.α_cabcar,WLSET.IΛ_PAR) : view(leaves[_i].BIO.α_cab,WLSET.IΛ_PAR);

        # convert energy to quantum unit for APAR and PPAR
        RADIATION._par_shaded  .= photon.(WLSET.Λ_PAR, view(RADIATION.e_sum_diffuse,WLSET.IΛ_PAR,_i)) .* 1000 ./ _tlai;
        RADIATION._par_sunlit  .= photon.(WLSET.Λ_PAR, view(RADIATION.e_sum_direct ,WLSET.IΛ_PAR,_i)) .* 1000 ./ _tlai;
        RADIATION._apar_shaded .= photon.(WLSET.Λ_PAR, view(RADIATION.e_net_diffuse,WLSET.IΛ_PAR,_i)) .* 1000 ./ _tlai;
        RADIATION._apar_sunlit .= photon.(WLSET.Λ_PAR, view(RADIATION.e_net_direct ,WLSET.IΛ_PAR,_i)) .* 1000 ./ _tlai;
        RADIATION._ppar_shaded .= RADIATION._apar_shaded .* _α_apar;
        RADIATION._ppar_sunlit .= RADIATION._apar_sunlit .* _α_apar;

        # PAR for leaves
        _Σ_par_dif = RADIATION._par_shaded' * WLSET.ΔΛ_PAR;
        _Σ_par_dir = RADIATION._par_sunlit' * WLSET.ΔΛ_PAR * _normi;
        RADIATION.par_shaded[_i] = _Σ_par_dif;
        RADIATION.par_sunlit[:,:,_i] .= OPTICS._abs_fs_fo .* _Σ_par_dir;
        RADIATION.par_sunlit[:,:,_i] .+= _Σ_par_dif;

        # APAR for leaves
        _Σ_apar_dif = RADIATION._apar_shaded' * WLSET.ΔΛ_PAR;
        _Σ_apar_dir = RADIATION._apar_sunlit' * WLSET.ΔΛ_PAR * _normi;
        RADIATION.apar_shaded[_i] = _Σ_apar_dif;
        RADIATION.apar_sunlit[:,:,_i] .= OPTICS._abs_fs_fo .* _Σ_apar_dir;
        RADIATION.apar_sunlit[:,:,_i] .+= _Σ_apar_dif;

        # PPAR for leaves
        _Σ_ppar_dif = RADIATION._apar_shaded' * WLSET.ΔΛ_PAR;
        _Σ_ppar_dir = RADIATION._apar_sunlit' * WLSET.ΔΛ_PAR * _normi;
        RADIATION.ppar_shaded[_i] = _Σ_ppar_dif;
        RADIATION.ppar_sunlit[:,:,_i] .= OPTICS._abs_fs_fo .* _Σ_ppar_dir;
        RADIATION.ppar_sunlit[:,:,_i] .+= _Σ_ppar_dif;
    end;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-10: migrate the function thermal_fluxes! from CanopyLayers
#     2022-Jun-10: update net lw radiation for leaves and soil
#
#######################################################################################################################################################################################################
"""

    canopy_radiation!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::FT, soil::Soil{FT}) where {FT<:AbstractFloat}

Updates canopy radiation profiles for longwave radiation, given
- `can` `HyperspectralMLCanopy` type struct
- `leaves` Vector of `Leaf`
- `rad` Incoming longwave radiation
- `soil` Bottom soil boundary layer
"""
canopy_radiation!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::FT, soil::Soil{FT}) where {FT<:AbstractFloat} = (
    @unpack N_LAYER, OPTICS, RADIATION = can;

    # 1. compute longwave radiation out from the leaves and soil
    for _i in eachindex(leaves)
        RADIATION.r_lw[_i] = K_STEFAN(FT) * OPTICS.ϵ[_i] * leaves[_i].t ^ 4;
    end;

    _r_lw_soil = K_STEFAN(FT) * (1 - soil.ρ_lw) * soil.t ^ 4;

    # 2. account for the longwave emission from bottom to up
    RADIATION._r_emit_up[end] = _r_lw_soil;

    for _i in N_LAYER:-1:1
        _r__ = OPTICS._ρ_lw[_i];
        _r_j = OPTICS.ρ_lw[_i+1];
        _t__ = OPTICS._τ_lw[_i];

        _dnorm = 1 - _r__ * _r_j;

        RADIATION._r_emit_down[_i] = (RADIATION._r_emit_up[_i+1] * _r__ + RADIATION.r_lw[_i]) / _dnorm;
        RADIATION._r_emit_up[_i] = RADIATION._r_emit_down[_i] * _r_j * _t__ + RADIATION._r_emit_up[_i+1] * _t__ + RADIATION.r_lw[_i];
    end;

    # 3. account for the longwave emission from up to bottom
    RADIATION.r_lw_down[1] = rad;

    for _i in 1:N_LAYER
        _r_i = OPTICS.ρ_lw[_i];
        _t_i = OPTICS.τ_lw[_i];

        RADIATION.r_lw_down[_i+1] = RADIATION.r_lw_down[_i] * _t_i + RADIATION._r_emit_down[_i];
        RADIATION.r_lw_up[_i+1] = RADIATION.r_lw_down[_i] * _r_i + RADIATION._r_emit_up[_i];
    end;

    RADIATION.r_lw_up[end] = RADIATION.r_lw_down[end] * soil.ρ_lw + _r_lw_soil;

    # 4. compute the net longwave radiation per canopy layer and soil
    for _i in 1:N_LAYER
        RADIATION.r_net_lw[_i] = (RADIATION.r_lw_down[_i] + RADIATION.r_lw_up[_i+1]) * (1 - OPTICS._ρ_lw[_i] - OPTICS._τ_lw[_i]) - 2* RADIATION.r_lw[_i];
    end;

    soil.r_net_lw = RADIATION.r_lw_down[end] * (1 - soil.ρ_lw) - _r_lw_soil;

    return nothing
);
