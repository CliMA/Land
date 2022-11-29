#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-09: migrate the function from CanopyLayers
#     2022-Jun-09: rename function to canopy_radiation!
#
#######################################################################################################################################################################################################
"""
This function updates canopy radiation profiles. The supported methods are to
- Update shortwave radiation profile for broadband or hyperspectral canopy
- Updates soil shortwave radiation profiles
- Update longwave radation profile for broadband or hyperspectral canopy
- Update radiation profile for SPAC

"""
function canopy_radiation! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-16: finalize the two leaf radiation
#     2022-Jun-21: make par and apar per leaf area
#     2022-Jun-21: redo the mean rates by weighing them by sunlit fraction
#     2022-Jun-25: finalize the soil shortwave reflection
#     2022-Jun-25: clean the calculations
#     2022-Jun-25: add method for longwave radation for two leaf RT
#     2022-Jun-29: use Leaves1D for the broadband RT
#     2022-Jun-29: use sunlit and shaded temperatures for the longwave out
#
#######################################################################################################################################################################################################
"""

    canopy_radiation!(can::BroadbandSLCanopy{FT}, leaf::Leaves1D{FT}, rad::BroadbandRadiation{FT}, soil::Soil{FT}) where {FT<:AbstractFloat}
    canopy_radiation!(can::BroadbandSLCanopy{FT}, leaf::Leaves1D{FT}, rad::FT, soil::Soil{FT}) where {FT<:AbstractFloat}

Updates shortwave or longwave radiation profiles, given
- `can` `HyperspectralMLCanopy` type struct
- `leaf` `Leaves1D` type struct
- `rad` Broadband shortwave or longwave radiation
- `soil` `Soil` type struct

"""
canopy_radiation!(can::BroadbandSLCanopy{FT}, leaf::Leaves1D{FT}, rad::BroadbandRadiation{FT}, soil::Soil{FT}) where {FT<:AbstractFloat} = (
    @unpack RADIATION = can;
    @unpack BIO = leaf;
    @unpack ALBEDO = soil;

    # compute the sunlit and shaded leaf area index of the entire canopy (adapted from Campbell 1998 equation 15.23, with clumping index)
    RADIATION.lai_sunlit = (1 - exp(-RADIATION.k_direct * can.lai * can.ci)) / RADIATION.k_direct;
    RADIATION.lai_shaded = can.lai - RADIATION.lai_sunlit;

    # theory for solar radiation
    #     q_sunlit(x) = q_sun_direct * k_direct + q_diffuse(x) + q_scatter(x) + sr(x)
    #     q_shaded(x) = q_diffuse(x) + q_scatter(x) + sr(x)
    #
    #     τ_directs(x) = exp(-sqrt(α) * k_direct * ci * x)
    #     τ_direct(x)  = exp(-k_direct * ci * x)
    #     τ_diffuse(x) = exp(-sqrt(α) * k_diffuse * ci * x)
    #
    #     q_directs(x) = τ_directs(x) * q_sun_direct   # include scattering, multiplied by sqrt(absorption)
    #     q_direct(x)  = τ_direct(x)  * q_sun_direct   # direct only
    #     q_diffuse(x) = τ_diffuse(x) * q_sun_diffuse
    #     q_scatter(x) = q_directs(x) - q_direct(x)
    #     q_soil_all   = q_directs(LAI) + q_diffuse(LAI)
    #     q_soil(x)    = τ_diffuse(LAI-x) * q_soil_all
    #
    # need to compute the leaf area weighted mean
    #     mean(q_direct)  = ∫q_direct(x)⋅dL / ∫dx
    #     mean(q_directs) = ∫q_directs(x)⋅dL / ∫dx
    #     mean(q_diffuse) = ∫q_diffuse(x)⋅dL / ∫dx
    #

    # when the radiation is from top to buttom
    @inline shaded_integral_top(k::FT, α::FT) where {FT<:AbstractFloat} = (
        return (1 - exp(-sqrt(α) * k * can.ci * can.lai)) / (sqrt(α) * k * can.ci) - (1 - exp(-(sqrt(α) * k + RADIATION.k_direct) * can.ci * can.lai)) / (sqrt(α) * k + RADIATION.k_direct)
    );
    @inline sunlit_integral_top(k::FT, α::FT) where {FT<:AbstractFloat} = (
        return (1 - exp(-(sqrt(α) * k + RADIATION.k_direct) * can.ci * can.lai)) / (sqrt(α) * k + RADIATION.k_direct)
    );

    # when the radiation is from bottom to top
    @inline shaded_integral_down(k::FT, α::FT) where {FT<:AbstractFloat} = (
        return exp(-sqrt(α) * k * can.ci * can.lai) / (RADIATION.k_direct + sqrt(α) * k) * (exp((RADIATION.k_direct + sqrt(α) * k) * can.ci * can.lai) - 1)
    );
    @inline sunlit_integral_down(k::FT, α::FT) where {FT<:AbstractFloat} = (
        return (1 - exp(-sqrt(α) * k * can.ci * can.lai)) / (sqrt(α) * k * can.ci) - shaded_integral_down(k, α)
    );

    # solar radiation reaching soil and the canopy (averaged)
    _soil_par     = rad.e_diffuse_par * exp(-sqrt(BIO.α_PAR) * RADIATION.k_diffuse * can.ci * can.lai) + rad.e_direct_par * exp(-sqrt(BIO.α_PAR) * RADIATION.k_direct * can.ci * can.lai);
    _soil_nir     = rad.e_diffuse_nir * exp(-sqrt(BIO.α_NIR) * RADIATION.k_diffuse * can.ci * can.lai) + rad.e_direct_nir * exp(-sqrt(BIO.α_NIR) * RADIATION.k_direct * can.ci * can.lai);
    _sunlit_s_par = _soil_par * (1 - ALBEDO.ρ_sw[1]) * sunlit_integral_down(RADIATION.k_diffuse, BIO.α_PAR);
    _shaded_s_par = _soil_par * (1 - ALBEDO.ρ_sw[1]) * shaded_integral_down(RADIATION.k_diffuse, BIO.α_PAR);
    _sunlit_s_nir = _soil_nir * (1 - ALBEDO.ρ_sw[2]) * sunlit_integral_down(RADIATION.k_diffuse, BIO.α_NIR);
    _shaded_s_nir = _soil_nir * (1 - ALBEDO.ρ_sw[2]) * shaded_integral_down(RADIATION.k_diffuse, BIO.α_NIR);

    # net radiation for soil
    ALBEDO.r_net_sw = _soil_par * (1 - ALBEDO.ρ_sw[1]) + _soil_nir * (1 - ALBEDO.ρ_sw[2]);

    # absorbed PAR for shaded leaves from solar radiation and soil reflection
    _shaded_q_diffuse     = rad.e_diffuse_par * shaded_integral_top(RADIATION.k_diffuse, BIO.α_PAR);
    _shaded_q_direct      = rad.e_direct_par  * shaded_integral_top(RADIATION.k_direct , FT(1));
    _shaded_q_directs     = rad.e_direct_par  * shaded_integral_top(RADIATION.k_direct , BIO.α_PAR);
    RADIATION.par_shaded  = (_shaded_q_diffuse + _shaded_q_directs - _shaded_q_direct + _shaded_s_par) / FT(0.235) / RADIATION.lai_shaded;
    RADIATION.apar_shaded = RADIATION.par_shaded * BIO.α_PAR;

    # absorbed PAR for sunlit leaves from solar radiation and soil reflection
    _sunlit_q_diffuse     = rad.e_diffuse_par * sunlit_integral_top(RADIATION.k_diffuse, BIO.α_PAR);
    _sunlit_q_direct      = rad.e_direct_par  * sunlit_integral_top(RADIATION.k_direct , FT(1));
    _sunlit_q_directs     = rad.e_direct_par  * sunlit_integral_top(RADIATION.k_direct , BIO.α_PAR);
    RADIATION.par_sunlit  = (rad.e_direct_par * RADIATION.k_direct + (_sunlit_q_diffuse + _sunlit_q_directs - _sunlit_q_direct + _sunlit_s_par) / FT(0.235) / RADIATION.lai_sunlit);
    RADIATION.apar_sunlit = RADIATION.par_sunlit * BIO.α_PAR;

    # absorbed NIR for shaded leaves
    _shaded_r_diffuse = rad.e_diffuse_nir * shaded_integral_top(RADIATION.k_diffuse, BIO.α_NIR);
    _shaded_r_direct  = rad.e_direct_nir  * shaded_integral_top(RADIATION.k_direct , FT(1));
    _shaded_r_directs = rad.e_direct_nir  * shaded_integral_top(RADIATION.k_direct , BIO.α_NIR);

    # absorbed NIR for sunlit leaves
    _sunlit_r_diffuse = rad.e_diffuse_nir * sunlit_integral_top(RADIATION.k_diffuse, BIO.α_NIR);
    _sunlit_r_direct  = rad.e_direct_nir  * sunlit_integral_top(RADIATION.k_direct , FT(1));
    _sunlit_r_directs = rad.e_direct_nir  * sunlit_integral_top(RADIATION.k_direct , BIO.α_NIR);

    # net radiation in-out for sunlit and shaded leaves for shortwave radation
    RADIATION.r_net_shaded = (_shaded_q_diffuse + _shaded_q_directs - _shaded_q_direct + _shaded_s_par) / RADIATION.lai_shaded * BIO.α_PAR +    # PAR region
                             (_shaded_r_diffuse + _shaded_r_directs - _shaded_r_direct + _shaded_s_nir) / RADIATION.lai_shaded * BIO.α_NIR;     # NIR region
    RADIATION.r_net_sunlit = (rad.e_direct_par * RADIATION.k_direct + (_sunlit_q_diffuse + _sunlit_q_directs - _sunlit_q_direct + _sunlit_s_par) / RADIATION.lai_sunlit) * BIO.α_PAR +
                             (rad.e_direct_nir * RADIATION.k_direct + (_sunlit_r_diffuse + _sunlit_r_directs - _sunlit_r_direct + _sunlit_s_nir) / RADIATION.lai_sunlit) * BIO.α_NIR;

    return nothing
);

canopy_radiation!(can::BroadbandSLCanopy{FT}, leaf::Leaves1D{FT}, rad::FT, soil::Soil{FT}) where {FT<:AbstractFloat} = (
    @unpack RADIATION = can;
    @unpack BIO = leaf;
    @unpack ALBEDO, LAYERS = soil;

    # theory for longwave radiation reaching soil
    #     soil_lw_in_solar = rad * τ_diffuse(LAI)
    #     soil_lw_in_sunlit(x) = ϵσT⁴(sunlit) * τ_direct(LAI-x)
    #     soil_lw_in_shaded(x) = ϵσT⁴(shaded) * τ_diffuse(LAI-x)
    #     soil_in_all = ∫(soil_lw_in_sunlit + soil_lw_in_shaded) + soil_lw_in_solar
    #
    # theory for longwave radiation reaching leaves
    #     leaf_lw_in_solar = rad * τ_diffuse(x)
    #     leaf_lw_in_soil = soil_out * τ_diffuse(LAI-x)

    # when the radiation is from top to buttom
    @inline shaded_integral_top(k::FT, α::FT) where {FT<:AbstractFloat} = (
        return (1 - exp(-sqrt(α) * k * can.ci * can.lai)) / (sqrt(α) * k * can.ci) - (1 - exp(-(sqrt(α) * k + RADIATION.k_direct) * can.ci * can.lai)) / (sqrt(α) * k + RADIATION.k_direct)
    );
    @inline sunlit_integral_top(k::FT, α::FT) where {FT<:AbstractFloat} = (
        return (1 - exp(-(sqrt(α) * k + RADIATION.k_direct) * can.ci * can.lai)) / (sqrt(α) * k + RADIATION.k_direct)
    );

    # when the radiation is from bottom to top
    @inline shaded_integral_down(k::FT, α::FT) where {FT<:AbstractFloat} = (
        return exp(-sqrt(α) * k * can.ci * can.lai) / (RADIATION.k_direct + sqrt(α) * k) * (exp((RADIATION.k_direct + sqrt(α) * k) * can.ci * can.lai) - 1)
    );
    @inline sunlit_integral_down(k::FT, α::FT) where {FT<:AbstractFloat} = (
        return (1 - exp(-sqrt(α) * k * can.ci * can.lai)) / (sqrt(α) * k * can.ci) - shaded_integral_down(k, α)
    );

    # longwave radiation reaching soil and the canopy (averaged)
    _soil_lw_in_solar  = rad * exp(-sqrt(BIO.ϵ_LW) * RADIATION.k_diffuse * can.ci * can.lai);
    _soil_lw_in_sunlit = BIO.ϵ_LW * K_STEFAN(FT) * leaf.t[1] ^ 4 * sunlit_integral_down(RADIATION.k_diffuse, BIO.ϵ_LW);
    _soil_lw_in_shaded = BIO.ϵ_LW * K_STEFAN(FT) * leaf.t[2] ^ 4 * shaded_integral_down(RADIATION.k_diffuse, BIO.ϵ_LW);
    _soil_lw_in        = _soil_lw_in_solar + _soil_lw_in_sunlit + _soil_lw_in_shaded;
    _soil_lw_out       = _soil_lw_in * ALBEDO.ρ_LW + (1 - ALBEDO.ρ_LW) * K_STEFAN(FT) + LAYERS[1].t ^ 4;
    ALBEDO.r_net_lw    = _soil_lw_in - _soil_lw_out;

    # lognwave radiation reaching sunlit and shaded leaves
    _shaded_in_solar  = rad * shaded_integral_top(RADIATION.k_diffuse, BIO.ϵ_LW);
    _sunlit_in_solar  = rad * sunlit_integral_top(RADIATION.k_diffuse, BIO.ϵ_LW);
    _shaded_in_soil   = _soil_lw_out * shaded_integral_down(RADIATION.k_diffuse, BIO.ϵ_LW);
    _sunlit_in_soil   = _soil_lw_out * sunlit_integral_down(RADIATION.k_diffuse, BIO.ϵ_LW);

    # lognwave radiation out from sunlit and shaded leaves
    _shaded_out_solar = BIO.ϵ_LW * K_STEFAN(FT) * leaf.t[1] ^ 4 * shaded_integral_top(RADIATION.k_diffuse, BIO.ϵ_LW);
    _sunlit_out_solar = BIO.ϵ_LW * K_STEFAN(FT) * leaf.t[2] ^ 4 * sunlit_integral_top(RADIATION.k_diffuse, BIO.ϵ_LW);

    # net radiation in-out for sunlit and shaded leaves for shortwave radation
    RADIATION.r_net_shaded += (_shaded_in_solar + _shaded_in_soil - _shaded_out_solar - _soil_lw_in_shaded) * BIO.ϵ_LW / RADIATION.lai_shaded;
    RADIATION.r_net_shaded += (_sunlit_in_solar + _sunlit_in_soil - _sunlit_out_solar - _soil_lw_in_sunlit) * BIO.ϵ_LW / RADIATION.lai_sunlit;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-14: make method work with broadband soil albedo struct
#     2022-Jun-14: allow method to use broadband PAR and NIR soil albedo values
#
#######################################################################################################################################################################################################
"""

    canopy_radiation!(can::HyperspectralMLCanopy{FT}, albedo::BroadbandSoilAlbedo{FT}) where {FT<:AbstractFloat}
    canopy_radiation!(can::HyperspectralMLCanopy{FT}, albedo::HyperspectralSoilAlbedo{FT}) where {FT<:AbstractFloat}

Updates soil shortwave radiation profiles, given
- `can` `HyperspectralMLCanopy` type struct
- `albedo` `BroadbandSoilAlbedo` or `HyperspectralSoilAlbedo` type soil albedo

"""
canopy_radiation!(can::HyperspectralMLCanopy{FT}, albedo::BroadbandSoilAlbedo{FT}) where {FT<:AbstractFloat} = (
    @unpack DIM_LAYER, OPTICS, RADIATION, WLSET = can;

    OPTICS._tmp_vec_λ[WLSET.IΛ_PAR] .= view(RADIATION.e_direct,WLSET.IΛ_PAR,DIM_LAYER+1) .* (1 .- albedo.ρ_sw[1]);
    OPTICS._tmp_vec_λ[WLSET.IΛ_NIR] .= view(RADIATION.e_direct,WLSET.IΛ_NIR,DIM_LAYER+1) .* (1 .- albedo.ρ_sw[2]);
    albedo.e_net_direct = OPTICS._tmp_vec_λ' * WLSET.ΔΛ / 1000;

    OPTICS._tmp_vec_λ[WLSET.IΛ_PAR] .= view(RADIATION.e_diffuse_down,WLSET.IΛ_PAR,DIM_LAYER+1) .* (1 .- albedo.ρ_sw[1]);
    OPTICS._tmp_vec_λ[WLSET.IΛ_NIR] .= view(RADIATION.e_diffuse_down,WLSET.IΛ_NIR,DIM_LAYER+1) .* (1 .- albedo.ρ_sw[2]);
    albedo.e_net_diffuse = OPTICS._tmp_vec_λ' * WLSET.ΔΛ / 1000;

    albedo.r_net_sw = albedo.e_net_direct + albedo.e_net_diffuse;

    return nothing
);

canopy_radiation!(can::HyperspectralMLCanopy{FT}, albedo::HyperspectralSoilAlbedo{FT}) where {FT<:AbstractFloat} = (
    @unpack DIM_LAYER, RADIATION, WLSET = can;

    albedo.e_net_direct .= view(RADIATION.e_direct,:,DIM_LAYER+1) .* (1 .- albedo.ρ_sw);
    albedo.e_net_diffuse .= view(RADIATION.e_diffuse_down,:,DIM_LAYER+1) .* (1 .- albedo.ρ_sw);
    albedo.r_net_sw = (albedo.e_net_direct' * WLSET.ΔΛ + albedo.e_net_diffuse' * WLSET.ΔΛ) / 1000;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-09: migrate the function from CanopyLayers
#     2022-Jun-10: rename PAR/APAR to APAR/PPAR to be more accurate
#     2022-Jun-10: add PAR calculation (before absorption)
#     2022-Jun-10: compute shortwave net radiation
#     2022-Jun-10: migrate the function thermal_fluxes! from CanopyLayers
#     2022-Jun-10: update net lw radiation for leaves and soil
#     2022-Jun-13: use DIM_LAYER instead of _end
#     2022-Jun-29: use Leaves2D for the hyperspectral RT
# Bug fix:
#     2022-Jul-15: sum by r_net_sw by the weights of sunlit and shaded fractions
#     2022-Jul-27: use _ρ_dd, _ρ_sd, _τ_dd, and _τ_sd for leaf energy absorption (typo when refactoring the code)
#     2022-Aug-30: fix par, apar, and ppar issues
#     2022-Nov-29: ppar for sunlit leaves should be scaled based on sunlit fraction
#
#######################################################################################################################################################################################################
"""

    canopy_radiation!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaves2D{FT}}, rad::HyperspectralRadiation{FT}, soil::Soil{FT}; APAR_CAR::Bool = true) where {FT<:AbstractFloat}
    canopy_radiation!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaves2D{FT}}, rad::FT, soil::Soil{FT}) where {FT<:AbstractFloat}

Updates canopy radiation profiles for shortwave or longwave radiation, given
- `can` `HyperspectralMLCanopy` type struct
- `leaves` Vector of `Leaves2D`
- `rad` Incoming shortwave or longwave radiation
- `soil` Bottom soil boundary layer
- `APAR_CAR` Whether carotenoid absorption is counted in PPAR, default is true

"""
canopy_radiation!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaves2D{FT}}, rad::HyperspectralRadiation{FT}, soil::Soil{FT}; APAR_CAR::Bool = true) where {FT<:AbstractFloat} = (
    @unpack DIM_LAYER, OPTICS, P_INCL, RADIATION, WLSET = can;
    @unpack ALBEDO = soil;
    _ilai = can.lai * can.ci / DIM_LAYER;
    _tlai = can.lai / DIM_LAYER;

    # 1. update upward and downward direct and diffuse radiation profiles
    RADIATION.e_direct[:,1] .= rad.e_direct;
    RADIATION.e_diffuse_down[:,1] .= rad.e_diffuse;

    for _i in 1:DIM_LAYER
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

    RADIATION.e_diffuse_up[:,end] .= view(OPTICS.ρ_sd,:,DIM_LAYER+1) .* view(RADIATION.e_direct,:,DIM_LAYER+1) .+ view(OPTICS.ρ_dd,:,DIM_LAYER+1) .* view(RADIATION.e_diffuse_down,:,DIM_LAYER+1);

    # 2. update the sunlit and shaded sum radiation and total absorbed radiation per layer and for soil
    for _i in 1:DIM_LAYER
        _a_s_i = view(RADIATION.e_net_direct  ,:,_i  );     # net absorbed direct radiation
        _a_d_i = view(RADIATION.e_net_diffuse ,:,_i  );     # net absorbed diffuse radiation
        _e_d_i = view(RADIATION.e_diffuse_down,:,_i  );     # downward diffuse radiation at upper boundary
        _e_s_i = view(RADIATION.e_direct      ,:,_i  );     # direct radiation at upper boundary
        _e_u_j = view(RADIATION.e_diffuse_up  ,:,_i+1);     # upward diffuse radiation at lower boundary
        _p_s_i = view(RADIATION.e_sum_direct  ,:,_i  );     # sum direct radiation
        _p_d_i = view(RADIATION.e_sum_diffuse ,:,_i  );     # sum diffuse radiation

        _r_dd__ = view(OPTICS._ρ_dd,:,_i);  # reflectance of the upper boundary (i)
        _r_sd__ = view(OPTICS._ρ_sd,:,_i);  # reflectance of the upper boundary (i)
        _t_dd__ = view(OPTICS._τ_dd,:,_i);  # transmittance of the layer (i)
        _t_sd__ = view(OPTICS._τ_sd,:,_i);  # transmittance of the layer (i)
        _t_ss__ = OPTICS._τ_ss;             # transmittance for directional->directional

        _p_s_i .= _e_s_i;
        _p_d_i .= _e_d_i .+ _e_u_j;

        _a_s_i .= _p_s_i .* (1 .- _t_ss__ .- _t_sd__ .- _r_sd__);
        _a_d_i .= _p_d_i .* (1 .- _t_dd__ .- _r_dd__);
    end;

    # 3. compute the spectra at the observer direction
    for _i in 1:DIM_LAYER
        _e_d_i = view(RADIATION.e_diffuse_down,:,_i);   # downward diffuse radiation at upper boundary
        _e_u_i = view(RADIATION.e_diffuse_up  ,:,_i);   # upward diffuse radiation at upper boundary
        _e_v_i = view(RADIATION.e_v,:,_i);

        _dob_i = view(OPTICS.σ_dob,:,_i);   # scattering coefficient backward for diffuse->observer
        _dof_i = view(OPTICS.σ_dof,:,_i);   # scattering coefficient forward for diffuse->observer
        _so__i = view(OPTICS.σ_so ,:,_i);   # bidirectional from solar to observer

        # _e_v_i .= (OPTICS.po[_i] .* _dob_i .* _e_d_i .+ OPTICS.po[_i] .* _dof_i .* _e_u_i .+ OPTICS.pso[_i] .* _so__i .* rad.e_direct) * _ilai;
        _e_v_i  .= OPTICS.po[_i] .* _dob_i .* _e_d_i;
        _e_v_i .+= OPTICS.po[_i] .* _dof_i .* _e_u_i;
        _e_v_i .+= OPTICS.pso[_i] .* _so__i .* rad.e_direct;
        _e_v_i .*= _ilai;
    end;
    RADIATION.e_v[:,end] .= OPTICS.po[end] .* view(RADIATION.e_diffuse_up,:,DIM_LAYER+1);

    for _i in eachindex(RADIATION.e_o)
        RADIATION.e_o[_i] = sum(view(RADIATION.e_v,_i,:)) / FT(pi);
    end;

    RADIATION.albedo .= RADIATION.e_o * FT(pi) ./ (rad.e_direct .+ rad.e_diffuse);

    # 4. compute net absorption for leaves and soil
    for _i in 1:DIM_LAYER
        _Σ_shaded = view(RADIATION.e_net_diffuse,:,_i)' * WLSET.ΔΛ / 1000 / _tlai;
        _Σ_sunlit = view(RADIATION.e_net_direct ,:,_i)' * WLSET.ΔΛ / 1000 / _tlai;
        RADIATION.r_net_sw_shaded[_i] = _Σ_shaded;
        RADIATION.r_net_sw_sunlit[_i] = _Σ_sunlit / OPTICS.p_sunlit[_i] + _Σ_shaded;
        RADIATION.r_net_sw[_i] = _Σ_shaded * (1 - OPTICS.p_sunlit[_i]) + _Σ_sunlit * OPTICS.p_sunlit[_i];
    end;

    canopy_radiation!(can, ALBEDO);

    # 5. compute top-of-canopy and leaf level PAR, APAR, and PPAR
    RADIATION._par_shaded .= photon.(WLSET.Λ_PAR, view(rad.e_diffuse,WLSET.IΛ_PAR)) .* 1000;
    RADIATION._par_sunlit .= photon.(WLSET.Λ_PAR, view(rad.e_direct ,WLSET.IΛ_PAR)) .* 1000;
    RADIATION.par_in_diffuse = RADIATION._par_shaded' * WLSET.ΔΛ_PAR;
    RADIATION.par_in_direct = RADIATION._par_sunlit' * WLSET.ΔΛ_PAR;
    RADIATION.par_in = RADIATION.par_in_diffuse + RADIATION.par_in_direct;

    mul!(OPTICS._tmp_vec_azi, OPTICS._abs_fs', P_INCL);
    _normi = 1 / mean(OPTICS._tmp_vec_azi);

    for _i in 1:DIM_LAYER
        _α_apar = APAR_CAR ? view(leaves[_i].BIO.α_cabcar,WLSET.IΛ_PAR) : view(leaves[_i].BIO.α_cab,WLSET.IΛ_PAR);

        # convert energy to quantum unit for APAR and PPAR
        RADIATION._par_shaded  .= photon.(WLSET.Λ_PAR, view(RADIATION.e_sum_diffuse,WLSET.IΛ_PAR,_i)) .* 1000;
        RADIATION._par_sunlit  .= photon.(WLSET.Λ_PAR, view(RADIATION.e_sum_direct ,WLSET.IΛ_PAR,_i)) .* 1000;
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
        _Σ_ppar_dif = RADIATION._ppar_shaded' * WLSET.ΔΛ_PAR;
        _Σ_ppar_dir = RADIATION._ppar_sunlit' * WLSET.ΔΛ_PAR * _normi;
        leaves[_i].ppar_shaded  = _Σ_ppar_dif;
        leaves[_i].ppar_sunlit .= OPTICS._abs_fs_fo .* _Σ_ppar_dir ./ OPTICS.p_sunlit[_i] .+ _Σ_ppar_dif;
    end;

    return nothing
);

canopy_radiation!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaves2D{FT}}, rad::FT, soil::Soil{FT}) where {FT<:AbstractFloat} = (
    @unpack DIM_LAYER, OPTICS, RADIATION = can;
    @unpack ALBEDO, LAYERS = soil;

    # 1. compute longwave radiation out from the leaves and soil
    for _i in eachindex(leaves)
        RADIATION.r_lw[_i] = K_STEFAN(FT) * OPTICS.ϵ[_i] * leaves[_i].t ^ 4;
    end;

    _r_lw_soil = K_STEFAN(FT) * (1 - ALBEDO.ρ_LW) * LAYERS[1].t ^ 4;

    # 2. account for the longwave emission from bottom to up
    RADIATION._r_emit_up[end] = _r_lw_soil;

    for _i in DIM_LAYER:-1:1
        _r__ = OPTICS._ρ_lw[_i];
        _r_j = OPTICS.ρ_lw[_i+1];
        _t__ = OPTICS._τ_lw[_i];

        _dnorm = 1 - _r__ * _r_j;

        RADIATION._r_emit_down[_i] = (RADIATION._r_emit_up[_i+1] * _r__ + RADIATION.r_lw[_i]) / _dnorm;
        RADIATION._r_emit_up[_i] = RADIATION._r_emit_down[_i] * _r_j * _t__ + RADIATION._r_emit_up[_i+1] * _t__ + RADIATION.r_lw[_i];
    end;

    # 3. account for the longwave emission from up to bottom
    RADIATION.r_lw_down[1] = rad;

    for _i in 1:DIM_LAYER
        _r_i = OPTICS.ρ_lw[_i];
        _t_i = OPTICS.τ_lw[_i];

        RADIATION.r_lw_down[_i+1] = RADIATION.r_lw_down[_i] * _t_i + RADIATION._r_emit_down[_i];
        RADIATION.r_lw_up[_i+1] = RADIATION.r_lw_down[_i] * _r_i + RADIATION._r_emit_up[_i];
    end;

    RADIATION.r_lw_up[end] = RADIATION.r_lw_down[end] * ALBEDO.ρ_LW + _r_lw_soil;

    # 4. compute the net longwave radiation per canopy layer and soil
    for _i in 1:DIM_LAYER
        RADIATION.r_net_lw[_i] = (RADIATION.r_lw_down[_i] + RADIATION.r_lw_up[_i+1]) * (1 - OPTICS._ρ_lw[_i] - OPTICS._τ_lw[_i]) - 2 * RADIATION.r_lw[_i];
    end;

    ALBEDO.r_net_lw = RADIATION.r_lw_down[end] * (1 - ALBEDO.ρ_LW) - _r_lw_soil;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-29: add method for SPAC
#     2022-Jul-28: update soil albedo at the very first step
#
#######################################################################################################################################################################################################
"""

    canopy_radiation!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat}

Updates canopy radiation profiles for shortwave and longwave radiation, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, `MonoMLTreeSPAC` type SPAC

"""
canopy_radiation!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat} = (
    @unpack ANGLES, CANOPY, LEAVES, RAD_LW, RAD_SW, SOIL = spac;

    soil_albedo!(CANOPY, SOIL);
    canopy_optical_properties!(CANOPY, ANGLES);
    canopy_optical_properties!(CANOPY, LEAVES, SOIL);
    canopy_radiation!(CANOPY, LEAVES, RAD_SW, SOIL; APAR_CAR = LEAVES[1].APAR_CAR);
    canopy_radiation!(CANOPY, LEAVES, RAD_LW, SOIL);

    return nothing
);
