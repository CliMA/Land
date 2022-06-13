#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-10: migrate the function from CanopyLayers
#     2022-Jun-10: rename function to canopy_fluorescence!
#
#######################################################################################################################################################################################################
"""
This function updates canopy fluorescence profiles. The supported methods include

$(METHODLIST)

"""
function canopy_fluorescence! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-10: migrate the function from CanopyLayers
#
#######################################################################################################################################################################################################
"""

    canopy_fluorescence!(can::HyperspectralMLCanopy{FT}; ϕ_photon::Bool = true) where {FT<:AbstractFloat}

Updates canopy radiation profiles for shortwave radiation, given
- `can` `HyperspectralMLCanopy` type struct
- `ϕ_photon` If true (default), convert photon to photon when computing SIF; otherwise, convert energy to energy
"""
canopy_fluorescence!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::HyperspectralRadiation{FT}; ϕ_photon::Bool = true) where {FT<:AbstractFloat} = (
    @unpack N_LAYER, OPTICS, P_INCL, RADIATION, WLS = can;
    _ilai = can.lai * can.ci / N_LAYER;

    # function tp weight matrices by inclination angles
    @inline lidf_weight(mat_0, mat_1) = (
        OPTICS._tmp_mat_incl_azi_1 .= mat_0 .* mat_1;
        mul!(OPTICS._tmp_vec_azi, OPTICS._tmp_mat_incl_azi_1', P_INCL);

        return mean(OPTICS._tmp_vec_azi)
    );

    # 1. compute SIF emissions for different layers
    for _i in 1:N_LAYER
        OPTICS._mat⁺ .= (leaves[_i].BIO.mat_b .+ leaves[_i].BIO.mat_f) ./ 2;
        OPTICS._mat⁻ .= (leaves[_i].BIO.mat_b .- leaves[_i].BIO.mat_f) ./ 2;

        # integrate the energy in each wave length bins
        OPTICS._tmp_vec_sife_1 .= view(RADIATION.e_direct      ,WLS.IΛ_SIFE,_i) .* WLS.ΔΛ_SIFE;
        OPTICS._tmp_vec_sife_2 .= view(RADIATION.e_diffuse_down,WLS.IΛ_SIFE,_i) .* WLS.ΔΛ_SIFE;
        OPTICS._tmp_vec_sife_3 .= view(RADIATION.e_diffuse_up  ,WLS.IΛ_SIFE,_i) .* WLS.ΔΛ_SIFE;

        # convert the energy to photons
        OPTICS._tmp_vec_sife_4 .= photon.(WLS.Λ_SIFE, OPTICS._tmp_vec_sife_1);
        OPTICS._tmp_vec_sife_5 .= photon.(WLS.Λ_SIFE, OPTICS._tmp_vec_sife_2);
        OPTICS._tmp_vec_sife_6 .= photon.(WLS.Λ_SIFE, OPTICS._tmp_vec_sife_3);

        # determine which ones to use depending on ϕ_photon
        if ϕ_photon
            _e_dir, _e_dif_down, _e_dif_up = OPTICS._tmp_vec_sife_4, OPTICS._tmp_vec_sife_5, OPTICS._tmp_vec_sife_6;
        else
            _e_dir, _e_dif_down, _e_dif_up = OPTICS._tmp_vec_sife_1, OPTICS._tmp_vec_sife_2, OPTICS._tmp_vec_sife_3;
        end;

        # convert the excitation radiation to fluorescence components
        mul!(OPTICS._tmp_vec_sif_1, OPTICS._mat⁺, _e_dir);          # SIF component from direct light (before scaling)
        mul!(OPTICS._tmp_vec_sif_2, OPTICS._mat⁻, _e_dir);          # SIF component from direct light (before scaling)
        mul!(OPTICS._tmp_vec_sif_3, OPTICS._mat⁺, _e_dif_down);     # SIF component from downward diffuse light for backward (before scaling)
        mul!(OPTICS._tmp_vec_sif_4, OPTICS._mat⁻, _e_dif_down);     # SIF component from downward diffuse light for backward (before scaling)
        mul!(OPTICS._tmp_vec_sif_5, OPTICS._mat⁺, _e_dif_up);       # SIF component from upward diffuse light for forward (before scaling)
        mul!(OPTICS._tmp_vec_sif_6, OPTICS._mat⁻, _e_dif_up);       # SIF component from upward diffuse light for forward (before scaling)

        # add up the fluorescence at various wavelength bins for sunlit and (up- and down-ward) diffuse SIF
        _ϕ_sunlit = view(RADIATION.ϕ_sunlit,:,:,_i);
        _ϕ_shaded = RADIATION.ϕ_shaded[_i];

        # upward and downward SIF from direct and diffuse radiation per leaf area
        RADIATION._s_shaded_up   .= OPTICS._tmp_vec_sif_3 .* lidf_weight(_ϕ_shaded, 1)              .+ OPTICS._tmp_vec_sif_4 .* lidf_weight(_ϕ_shaded, can._COS²_Θ_INCL_AZI) .+
                                    OPTICS._tmp_vec_sif_5 .* lidf_weight(_ϕ_shaded, 1)              .- OPTICS._tmp_vec_sif_6 .* lidf_weight(_ϕ_shaded, can._COS²_Θ_INCL_AZI);
        RADIATION._s_shaded_down .= OPTICS._tmp_vec_sif_3 .* lidf_weight(_ϕ_shaded, 1)              .- OPTICS._tmp_vec_sif_4 .* lidf_weight(_ϕ_shaded, can._COS²_Θ_INCL_AZI) .+
                                    OPTICS._tmp_vec_sif_5 .* lidf_weight(_ϕ_shaded, 1)              .+ OPTICS._tmp_vec_sif_6 .* lidf_weight(_ϕ_shaded, can._COS²_Θ_INCL_AZI);
        RADIATION._s_sunlit_up   .= OPTICS._tmp_vec_sif_1 .* lidf_weight(_ϕ_sunlit, OPTICS._abs_fs) .+ OPTICS._tmp_vec_sif_2 .* lidf_weight(_ϕ_sunlit, OPTICS._fs_cos_θ_incl) .+
                                    OPTICS._tmp_vec_sif_3 .* lidf_weight(_ϕ_sunlit, 1)              .+ OPTICS._tmp_vec_sif_4 .* lidf_weight(_ϕ_sunlit, can._COS²_Θ_INCL_AZI) .+
                                    OPTICS._tmp_vec_sif_5 .* lidf_weight(_ϕ_sunlit, 1)              .- OPTICS._tmp_vec_sif_6 .* lidf_weight(_ϕ_sunlit, can._COS²_Θ_INCL_AZI);
        RADIATION._s_sunlit_down .= OPTICS._tmp_vec_sif_1 .* lidf_weight(_ϕ_sunlit, OPTICS._abs_fs) .- OPTICS._tmp_vec_sif_2 .* lidf_weight(_ϕ_sunlit, OPTICS._fs_cos_θ_incl) .+
                                    OPTICS._tmp_vec_sif_3 .* lidf_weight(_ϕ_sunlit, 1)              .- OPTICS._tmp_vec_sif_4 .* lidf_weight(_ϕ_sunlit, can._COS²_Θ_INCL_AZI) .+
                                    OPTICS._tmp_vec_sif_5 .* lidf_weight(_ϕ_sunlit, 1)              .+ OPTICS._tmp_vec_sif_6 .* lidf_weight(_ϕ_sunlit, can._COS²_Θ_INCL_AZI);

        # total emitted SIF for upward and downward direction
        RADIATION.s_layer_down[:,_i] .= _ilai .* OPTICS.ps[_i] .* RADIATION._s_sunlit_down .+ _ilai .* (1 - OPTICS.ps[_i]) .* RADIATION._s_shaded_down;
        RADIATION.s_layer_up[:,_i]   .= _ilai .* OPTICS.ps[_i] .* RADIATION._s_sunlit_up   .+ _ilai .* (1 - OPTICS.ps[_i]) .* RADIATION._s_shaded_up;
    end;

    # 2. account for the SIF emission from bottom to up
    RADIATION._s_emit_up[:,end] .= 0;

    for _i in N_LAYER:-1:1
        _r__ = view(OPTICS._ρ_dd,WLS.IΛ_SIF,_i  );  # reflectance without correction
        _r_j = view(OPTICS.ρ_dd ,WLS.IΛ_SIF,_i+1);  # reflectance of the upper boundary (i) for SIF
        _t_i = view(OPTICS.τ_dd ,WLS.IΛ_SIF,_i  );  # transmittance of the layer (i) for SIF

        _s_d_i = view(RADIATION._s_emit_down,:,_i  );   # downward SIF from the layer
        _s_u_i = view(RADIATION._s_emit_up  ,:,_i  );   # upward SIF from the layer
        _s_u_j = view(RADIATION._s_emit_up  ,:,_i+1);   # upward SIF from the lower layer
        _f_d_i = view(RADIATION.s_layer_down,:,_i  );   # downward emitted SIF from layer i
        _f_u_i = view(RADIATION.s_layer_up  ,:,_i  );   # downward emitted SIF from layer i

        OPTICS._tmp_vec_sif_1 .= 1 .- _r__ .* _r_j;

        _s_d_i .= (_f_d_i .+ _s_u_j .* _r__) ./ OPTICS._tmp_vec_sif_1;
        _s_u_i .= _f_u_i .+ _s_u_j .* _t_i .+ _s_d_i .* _r_j .* _t_i;
    end;

    # 3. account for the SIF emission from up to bottom
    RADIATION.sif_down[:,1] .= 0;

    for _i in 1:N_LAYER
        _r_i = view(OPTICS.ρ_dd,WLS.IΛ_SIF,_i);     # reflectance of the layer (i) for SIF
        _t_i = view(OPTICS.τ_dd,WLS.IΛ_SIF,_i);     # transmittance of the layer (i) for SIF

        _s_d_i = view(RADIATION._s_emit_down,:,_i  );   # downward SIF from the layer
        _s_u_i = view(RADIATION._s_emit_up  ,:,_i  );   # upward SIF from the layer

        _a_d_i = view(RADIATION.sif_down,:,_i  );
        _a_d_j = view(RADIATION.sif_down,:,_i+1);
        _a_u_i = view(RADIATION.sif_up  ,:,_i  );

        _a_d_j .= _a_d_i .* _t_i .+ _s_d_i;
        _a_u_i .= _a_d_i .* _r_i .+ _s_u_i;
    end;

    _end = lastindex(OPTICS.ρ_dd, 2);
    RADIATION.sif_up[:,end] .= view(RADIATION.sif_down,:,_end) .* view(OPTICS.ρ_dd,WLS.IΛ_SIF,_end) .+ view(RADIATION._s_emit_up,:,_end);

    return nothing
);
