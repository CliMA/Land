
function canopy_radiation! end



canopy_radiation!(can::HyperspectralMLCanopy{FT}, rad::HyperspectralRadiation{FT}, soil::Soil{FT}) where {FT<:AbstractFloat} = (
    @unpack N_LAYER, OPTICS, RADIATION = can;
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

    # update the sunlit and shaded total radiation per layer and for soil
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
        _e_d_i = view(RADIATION.e_diffuse_down,:,_i  );     # downward diffuse radiation at upper boundary
        _e_u_i = view(RADIATION.e_diffuse_up  ,:,_i  );     # upward diffuse radiation at upper boundary

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

    # compute soil net absorption

    return nothing
);
