
function canopy_radiation! end



canopy_radiation!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::HyperspectralRadiation{FT}, soil::Soil{FT}) where {FT<:AbstractFloat} = (
    @unpack APAR_CAR, N_LAYER, OPTICS, P_INCL, RADIATION, WLSET = can;
    _ilai = can.lai * can.ci / N_LAYER;
    _tlai = can.lai / N_LAYER;

    # update upward and downward direct and diffuse radiation profiles
    RADIATION.e_direct[:,1] .= rad.e_direct;
    RADIATION.e_diffuse_down[:,1] .= rad.e_diffuse;

    for _i in 1:N_LAYER
        _e_d_i = view(RADIATION.e_diffuse_down,:,_i  );     # downward diffuse radiation at upper boundary
        _e_d_j = view(RADIATION.e_diffuse_down,:,_i+1);     # downward diffuse radiation at lower boundary
        _e_s_i = view(RADIATION.e_direct      ,:,_i  );     # direct radiation at upper boundary
        _e_s_j = view(RADIATION.e_direct      ,:,_i+1);     # direct radiation at lower boundary
        _e_u_i = view(RADIATION.e_diffuse_up  ,:,_i  );     # upward diffuse radiation at upper boundary

        _r_dd_i = view(OPTICS.ρ_dd ,:,_i);  # reflectance of the upper boundary (i)
        _r_sd_i = view(OPTICS.ρ_sd ,:,_i);  # reflectance of the upper boundary (i)
        _t_dd_i = view(OPTICS.τ_dd ,:,_i);  # transmittance of the layer (i)
        _t_sd_i = view(OPTICS.τ_sd ,:,_i);  # transmittance of the layer (i)
        _t_ss__ = OPTICS._τ_ss;             # transmittance for directional->directional

        _e_s_j .= _t_ss__ .* _e_s_i;
        _e_d_j .= _t_sd_i .* _e_s_i .+ _t_dd_i .* _e_d_i;
        _e_u_i .= _r_sd_i .* _e_s_i .+ _r_dd_i .* _e_d_i;
    end;

    _end = lastindex(OPTICS.ρ_sd, 2);
    RADIATION.e_diffuse_up[:,end] = view(OPTICS.ρ_sd,:,_end) .* view(RADIATION.e_direct,:,_end) .+ view(OPTICS.ρ_dd,:,_end) .* view(RADIATION.e_diffuse_down,:,_end);

    # update the sunlit and shaded total absorbed radiation per layer and for soil
    for _i in 1:N_LAYER
        _a_s_i = view(RADIATION.e_net_direct  ,:,_i  );     # net absorbed direct radiation
        _a_d_i = view(RADIATION.e_net_diffuse ,:,_i  );     # net absorbed diffuse radiation
        _e_d_i = view(RADIATION.e_diffuse_down,:,_i  );     # downward diffuse radiation at upper boundary
        _e_s_i = view(RADIATION.e_direct      ,:,_i  );     # direct radiation at upper boundary
        _e_u_j = view(RADIATION.e_diffuse_up  ,:,_i+1);     # upward diffuse radiation at lower boundary

        _r_dd_i = view(OPTICS.ρ_dd ,:,_i);  # reflectance of the upper boundary (i)
        _r_sd_i = view(OPTICS.ρ_sd ,:,_i);  # reflectance of the upper boundary (i)
        _t_dd_i = view(OPTICS.τ_dd ,:,_i);  # transmittance of the layer (i)
        _t_sd_i = view(OPTICS.τ_sd ,:,_i);  # transmittance of the layer (i)
        _t_ss__ = OPTICS._τ_ss;             # transmittance for directional->directional

        _a_s_i .= _e_s_i .* (1 .- _t_ss__ .- _t_sd_i .- _r_sd_i);
        _a_d_i .= (_e_d_i .+ _e_u_j) .* (1 .- _t_dd_i .- _r_dd_i);
    end;

    soil.e_net_direct .= view(RADIATION.e_direct,:,_end) .* (1 .- soil.ρ_SW);
    soil.e_net_diffuse .= view(RADIATION.e_diffuse_down,:,_end) .* (1 .- soil.ρ_SW);

    # compute the spectra at the observer direction
    for _i in 1:N_LAYER
        _e_d_i = view(RADIATION.e_diffuse_down,:,_i);   # downward diffuse radiation at upper boundary
        _e_u_i = view(RADIATION.e_diffuse_up  ,:,_i);   # upward diffuse radiation at upper boundary

        _dob_i = view(OPTICS.σ_dob,:,_i);   # scattering coefficient backward for diffuse->observer
        _dof_i = view(OPTICS.σ_dob,:,_i);   # scattering coefficient forward for diffuse->observer
        _so__i = view(OPTICS.σ_so ,:,_i);   # bidirectional from solar to observer

        RADIATION.e_v[:,_i] .= (OPTICS.po[_i] .* _dob_i .* _e_d_i .+ OPTICS.po[_i] .* _dof_i .* _e_u_i .+ OPTICS.poo[_i] .* _so__i .* rad.e_direct) * _ilai;
    end;
    RADIATION.e_v[:,end] .= OPTICS.po[end] .* view(RADIATION.e_diffuse_up,:,_end);

    for _i in eachindex(RADIATION.e_o)
        RADIATION.e_o[_i] = sum(view(RADIATION.e_o,_i,:)) / FT(pi);
    end;

    RADIATION.albedo .= RADIATION.e_o * FT(pi) ./ (rad.e_direct .+ rad.e_diffuse);

    # compute net absorption for soil
    soil.r_net = (soil.e_net_direct' * WLSET.ΔΛ + soil.e_net_diffuse' * WLSET.ΔΛ) / 1000;

    # compute leaf level PAR and APAR
    mul!(OPTICS._tmp_vec_azi, OPTICS._abs_fs', P_INCL);
    _normi = 1 / mean(OPTICS._tmp_vec_azi);

    for _i in 1:N_LAYER
        _α_apar = APAR_CAR ? view(leaves[_i].BIO.α_cabcar,WLSET.IΛ_PAR) : view(leaves[_i].BIO.α_cab,WLSET.IΛ_PAR);

        # convert energy to quantum unit for APAR and PPAR
        RADIATION._apar_shaded .= photon.(WLSET.Λ_PAR, view(RADIATION.e_net_diffuse, WLSET.IΛ_PAR, _i)) .* 1000 ./ _tlai;
        RADIATION._apar_sunlit .= photon.(WLSET.Λ_PAR, view(RADIATION.e_net_direct , WLSET.IΛ_PAR, _i)) .* 1000 ./ _tlai;
        RADIATION._ppar_shaded .= RADIATION._apar_shaded .* _α_apar;
        RADIATION._ppar_sunlit .= RADIATION._apar_sunlit .* _α_apar;

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
