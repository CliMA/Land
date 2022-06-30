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
#     2022-Jun-10: use more descriptive variable names
#     2022-Jun-13: finished migrating the SIF function
#     2022-Jun-13: fix documentation
#     2022-Jun-14: convert energy and photon back and forth if using photon mode
#     2022-Jun-29: use Leaves2D for the hyperspectral RT
#     2022-Jun-29: use ϕ_f in Leaves2D
#
#######################################################################################################################################################################################################
"""

    canopy_fluorescence!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaves2D{FT}}; ϕ_photon::Bool = true) where {FT<:AbstractFloat}

Updates canopy radiation profiles for shortwave radiation, given
- `can` `HyperspectralMLCanopy` type struct
- `leaves` Vector of `Leaves2D`
- `ϕ_photon` If true (default), convert photon to photon when computing SIF; otherwise, convert energy to energy
"""
canopy_fluorescence!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaves2D{FT}}; ϕ_photon::Bool = true) where {FT<:AbstractFloat} = (
    @unpack N_LAYER, OPTICS, P_INCL, RADIATION, WLSET = can;
    _ilai = can.lai * can.ci / N_LAYER;

    # function to weight matrices by inclination angles
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
        OPTICS._tmp_vec_sife_1 .= view(RADIATION.e_direct      ,WLSET.IΛ_SIFE,1 ) .* WLSET.ΔΛ_SIFE;
        OPTICS._tmp_vec_sife_2 .= view(RADIATION.e_diffuse_down,WLSET.IΛ_SIFE,_i) .* WLSET.ΔΛ_SIFE;
        OPTICS._tmp_vec_sife_3 .= view(RADIATION.e_diffuse_up  ,WLSET.IΛ_SIFE,_i) .* WLSET.ΔΛ_SIFE;

        _e_dir, _e_dif_down, _e_dif_up = OPTICS._tmp_vec_sife_1, OPTICS._tmp_vec_sife_2, OPTICS._tmp_vec_sife_3;

        # determine which ones to use depending on ϕ_photon
        if ϕ_photon
            photon!(WLSET.Λ_SIFE, OPTICS._tmp_vec_sife_1);
            photon!(WLSET.Λ_SIFE, OPTICS._tmp_vec_sife_2);
            photon!(WLSET.Λ_SIFE, OPTICS._tmp_vec_sife_3);
        end;

        # convert the excitation radiation to fluorescence components
        mul!(OPTICS._tmp_vec_sif_1, OPTICS._mat⁺, _e_dir);          # SIF component from direct light (before scaling)
        mul!(OPTICS._tmp_vec_sif_2, OPTICS._mat⁻, _e_dir);          # SIF component from direct light (before scaling)
        mul!(OPTICS._tmp_vec_sif_3, OPTICS._mat⁺, _e_dif_down);     # SIF component from downward diffuse light for backward (before scaling)
        mul!(OPTICS._tmp_vec_sif_4, OPTICS._mat⁻, _e_dif_down);     # SIF component from downward diffuse light for backward (before scaling)
        mul!(OPTICS._tmp_vec_sif_5, OPTICS._mat⁺, _e_dif_up);       # SIF component from upward diffuse light for forward (before scaling)
        mul!(OPTICS._tmp_vec_sif_6, OPTICS._mat⁻, _e_dif_up);       # SIF component from upward diffuse light for forward (before scaling)

        # convert the SIF back to energy unit if ϕ_photon is true
        if ϕ_photon
            energy!(WLSET.Λ_SIF, OPTICS._tmp_vec_sif_1);
            energy!(WLSET.Λ_SIF, OPTICS._tmp_vec_sif_2);
            energy!(WLSET.Λ_SIF, OPTICS._tmp_vec_sif_3);
            energy!(WLSET.Λ_SIF, OPTICS._tmp_vec_sif_4);
            energy!(WLSET.Λ_SIF, OPTICS._tmp_vec_sif_5);
            energy!(WLSET.Λ_SIF, OPTICS._tmp_vec_sif_6);
        end;

        # add up the fluorescence at various wavelength bins for sunlit and (up- and down-ward) diffuse SIF
        _ϕ_sunlit = leaves[_i].ϕ_f_sunlit;
        _ϕ_shaded = leaves[_i].ϕ_f_shaded;

        # compute the weights
        _sh_1_ = lidf_weight(_ϕ_shaded, 1);
        _sl_1_ = lidf_weight(_ϕ_sunlit, 1);
        _sh_O_ = lidf_weight(_ϕ_shaded, OPTICS._abs_fo);
        _sl_O_ = lidf_weight(_ϕ_sunlit, OPTICS._abs_fo);
        _sl_S_ = lidf_weight(_ϕ_sunlit, OPTICS._abs_fs);
        _sh_oθ = lidf_weight(_ϕ_shaded, OPTICS._fo_cos_θ_incl);
        _sl_oθ = lidf_weight(_ϕ_sunlit, OPTICS._fo_cos_θ_incl);
        _sl_sθ = lidf_weight(_ϕ_sunlit, OPTICS._fs_cos_θ_incl);
        _sl_SO = lidf_weight(_ϕ_sunlit, OPTICS._abs_fs_fo);
        _sl_so = lidf_weight(_ϕ_sunlit, OPTICS._fs_fo);
        _sh_θ² = lidf_weight(_ϕ_shaded, can._COS²_Θ_INCL_AZI);
        _sl_θ² = lidf_weight(_ϕ_sunlit, can._COS²_Θ_INCL_AZI);

        # upward and downward SIF from direct and diffuse radiation per leaf area
        RADIATION._s_shaded_up   .= OPTICS._tmp_vec_sif_3 .* _sh_1_ .+ OPTICS._tmp_vec_sif_4 .* _sh_θ² .+
                                    OPTICS._tmp_vec_sif_5 .* _sh_1_ .- OPTICS._tmp_vec_sif_6 .* _sh_θ²;
        RADIATION._s_shaded_down .= OPTICS._tmp_vec_sif_3 .* _sh_1_ .- OPTICS._tmp_vec_sif_4 .* _sh_θ² .+
                                    OPTICS._tmp_vec_sif_5 .* _sh_1_ .+ OPTICS._tmp_vec_sif_6 .* _sh_θ²;
        RADIATION._s_sunlit_up   .= OPTICS._tmp_vec_sif_1 .* _sl_S_ .+ OPTICS._tmp_vec_sif_2 .* _sl_sθ .+
                                    OPTICS._tmp_vec_sif_3 .* _sl_1_ .+ OPTICS._tmp_vec_sif_4 .* _sl_θ² .+
                                    OPTICS._tmp_vec_sif_5 .* _sl_1_ .- OPTICS._tmp_vec_sif_6 .* _sl_θ²;
        RADIATION._s_sunlit_down .= OPTICS._tmp_vec_sif_1 .* _sl_S_ .- OPTICS._tmp_vec_sif_2 .* _sl_sθ .+
                                    OPTICS._tmp_vec_sif_3 .* _sl_1_ .- OPTICS._tmp_vec_sif_4 .* _sl_θ² .+
                                    OPTICS._tmp_vec_sif_5 .* _sl_1_ .+ OPTICS._tmp_vec_sif_6 .* _sl_θ²;

        # update the SIF cache for the observer direction (compute it here to save time)
        RADIATION._sif_obs_sunlit[:,_i] .= OPTICS._tmp_vec_sif_1 .* _sl_SO .+ OPTICS._tmp_vec_sif_2 .* _sl_so .+
                                           OPTICS._tmp_vec_sif_3 .* _sl_O_ .+ OPTICS._tmp_vec_sif_2 .* _sl_oθ .+
                                           OPTICS._tmp_vec_sif_5 .* _sl_O_ .- OPTICS._tmp_vec_sif_6 .* _sl_oθ;
        RADIATION._sif_obs_shaded[:,_i] .= OPTICS._tmp_vec_sif_3 .* _sh_O_ .+ OPTICS._tmp_vec_sif_2 .* _sh_oθ .+
                                           OPTICS._tmp_vec_sif_5 .* _sh_O_ .- OPTICS._tmp_vec_sif_6 .* _sh_oθ;

        # total emitted SIF for upward and downward direction
        RADIATION.s_layer_down[:,_i] .= _ilai .* OPTICS.p_sunlit[_i] .* RADIATION._s_sunlit_down .+ _ilai .* (1 - OPTICS.p_sunlit[_i]) .* RADIATION._s_shaded_down;
        RADIATION.s_layer_up[:,_i]   .= _ilai .* OPTICS.p_sunlit[_i] .* RADIATION._s_sunlit_up   .+ _ilai .* (1 - OPTICS.p_sunlit[_i]) .* RADIATION._s_shaded_up;
    end;

    # 2. account for the SIF emission from bottom to up
    RADIATION._s_emit_up[:,end] .= 0;

    for _i in N_LAYER:-1:1
        _r__ = view(OPTICS._ρ_dd,WLSET.IΛ_SIF,_i  );  # reflectance without correction
        _r_j = view(OPTICS.ρ_dd ,WLSET.IΛ_SIF,_i+1);  # reflectance of the upper boundary (i) for SIF
        _t_i = view(OPTICS.τ_dd ,WLSET.IΛ_SIF,_i  );  # transmittance of the layer (i) for SIF

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
        _r_i = view(OPTICS.ρ_dd,WLSET.IΛ_SIF,_i);   # reflectance of the layer (i) for SIF
        _t_i = view(OPTICS.τ_dd,WLSET.IΛ_SIF,_i);   # transmittance of the layer (i) for SIF

        _s_d_i = view(RADIATION._s_emit_down,:,_i);     # downward SIF from the layer
        _s_u_i = view(RADIATION._s_emit_up  ,:,_i);     # upward SIF from the layer

        _a_d_i = view(RADIATION.sif_down,:,_i  );
        _a_d_j = view(RADIATION.sif_down,:,_i+1);
        _a_u_i = view(RADIATION.sif_up  ,:,_i  );

        _a_d_j .= _a_d_i .* _t_i .+ _s_d_i;
        _a_u_i .= _a_d_i .* _r_i .+ _s_u_i;
    end;

    RADIATION.sif_up[:,end] .= view(RADIATION.sif_down,:,N_LAYER+1) .* view(OPTICS.ρ_dd,WLSET.IΛ_SIF,N_LAYER+1) .+ view(RADIATION._s_emit_up,:,N_LAYER+1);

    # 4. compute SIF from the observer direction
    OPTICS._tmp_vec_layer .= (view(OPTICS.pso,1:N_LAYER) .+ view(OPTICS.pso,2:N_LAYER+1)) ./ 2 .* _ilai ./ FT(pi);
    mul!(RADIATION.sif_obs_sunlit, RADIATION._sif_obs_sunlit, OPTICS._tmp_vec_layer);

    OPTICS._tmp_vec_layer .= (view(OPTICS.po,1:N_LAYER) .+ view(OPTICS.po,2:N_LAYER+1) .- view(OPTICS.pso,1:N_LAYER) .- view(OPTICS.pso,2:N_LAYER+1)) ./ 2 .* _ilai ./ FT(pi);
    mul!(RADIATION.sif_obs_shaded, RADIATION._sif_obs_shaded, OPTICS._tmp_vec_layer);

    RADIATION._sif_obs_scatter .= view(OPTICS.σ_ddb,WLSET.IΛ_SIF,:) .* view(RADIATION.sif_down,:,1:N_LAYER) .+ view(OPTICS.σ_ddf,WLSET.IΛ_SIF,:) .* view(RADIATION.sif_up,:,1:N_LAYER);
    OPTICS._tmp_vec_layer .= (view(OPTICS.po,1:N_LAYER) .+ view(OPTICS.po,2:N_LAYER+1)) ./ 2 .* _ilai ./ FT(pi);
    mul!(RADIATION.sif_obs_scatter, RADIATION._sif_obs_scatter, OPTICS._tmp_vec_layer);

    RADIATION.sif_obs .= RADIATION.sif_obs_sunlit .+ RADIATION.sif_obs_shaded .+ RADIATION.sif_obs_scatter .+ view(RADIATION.sif_up,:,N_LAYER+1) ./ OPTICS.po[end] ./ FT(pi);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-29: add method for SPAC
#
#######################################################################################################################################################################################################
"""

    canopy_fluorescence!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat}

Updates canopy fluorescence, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, `MonoMLTreeSPAC` type SPAC
"""
canopy_fluorescence!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat} = (
    @unpack CANOPY, LEAVES, Φ_PHOTON = spac;

    canopy_fluorescence!(CANOPY, LEAVES; ϕ_photon = Φ_PHOTON);

    return nothing
);
