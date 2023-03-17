#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-07: migrate the function from CanopyLayers
#     2022-Jun-09: rename function to canopy_optical_properties!
#
#######################################################################################################################################################################################################
"""
This function updates canopy optical properties for canopy. The supported methods are to
- Update the extinction coefficients
- Update the soil boundary conditions (not public function)
- Update scattering coefficient matrices

"""
function canopy_optical_properties! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-07: migrate the function from CanopyLayers
#     2022-Jun-09: rename function to canopy_optical_properties!
#     2022-Jun-09: update p_sunlit
#     2023-Mar-11: add code to account for the case of LAI == 0
#
#######################################################################################################################################################################################################
"""

    canopy_optical_properties!(can::HyperspectralMLCanopy{FT}, angles::SunSensorGeometry{FT}) where {FT<:AbstractFloat}

Updates canopy optical properties (extinction coefficients for direct and diffuse light) based on the SAIL model, given
- `can` `HyperspectralMLCanopy` type struct
- `angles` `SunSensorGeometry` type struct

"""
canopy_optical_properties!(can::HyperspectralMLCanopy{FT}, angles::SunSensorGeometry{FT}) where {FT<:AbstractFloat} = (
    (; DIM_LAYER, HOT_SPOT, OPTICS, P_INCL, Θ_AZI) = can;

    if can.lai == 0
        return nothing
    end;

    # 1. update the canopy optical properties related to extinction and scattering coefficients
    extinction_scattering_coefficients!(can, angles);

    OPTICS.ko  = P_INCL' * OPTICS._ko;
    OPTICS.ks  = P_INCL' * OPTICS._ks;
    OPTICS.sob = P_INCL' * OPTICS._sb;
    OPTICS.sof = P_INCL' * OPTICS._sf;
    OPTICS._bf = P_INCL' * can._COS²_Θ_INCL;

    OPTICS.sdb = (OPTICS.ks + OPTICS._bf) / 2;
    OPTICS.sdf = (OPTICS.ks - OPTICS._bf) / 2;
    OPTICS.dob = (OPTICS.ko + OPTICS._bf) / 2;
    OPTICS.dof = (OPTICS.ko - OPTICS._bf) / 2;
    OPTICS.ddb = (1 + OPTICS._bf) / 2;
    OPTICS.ddf = (1 - OPTICS._bf) / 2;

    # 2. update the matrices fs and fo
    OPTICS._cos_θ_azi_raa .= cosd.(Θ_AZI .- (angles.vaa - angles.saa));
    mul!(OPTICS._tmp_mat_incl_azi_1, OPTICS._Co, can._1_AZI');
    mul!(OPTICS._tmp_mat_incl_azi_2, OPTICS._So, OPTICS._cos_θ_azi_raa');
    OPTICS.fo .= (OPTICS._tmp_mat_incl_azi_1 .+ OPTICS._tmp_mat_incl_azi_2) ./ cosd(angles.vza);
    OPTICS._abs_fo .= abs.(OPTICS.fo);

    mul!(OPTICS._tmp_mat_incl_azi_1, OPTICS._Cs, can._1_AZI');
    mul!(OPTICS._tmp_mat_incl_azi_2, OPTICS._Ss, can._COS_Θ_AZI');
    OPTICS.fs .= (OPTICS._tmp_mat_incl_azi_1 .+ OPTICS._tmp_mat_incl_azi_2) ./ cosd(angles.sza);
    OPTICS._abs_fs .= abs.(OPTICS.fs);

    OPTICS._fo_cos_θ_incl .= OPTICS.fo .* can._COS²_Θ_INCL_AZI;
    OPTICS._fs_cos_θ_incl .= OPTICS.fs .* can._COS²_Θ_INCL_AZI;
    OPTICS._fs_fo .= OPTICS.fs .* OPTICS.fo;
    OPTICS._abs_fs_fo .= abs.(OPTICS._fs_fo);

    # 3. update the viewing fraction ps, po, pso, and p_sunlit
    _fac_s = (1 - exp(-OPTICS.ks * can.ci * can.lai / DIM_LAYER)) / (OPTICS.ks * can.ci * can.lai / DIM_LAYER);
    _fac_o = (1 - exp(-OPTICS.ko * can.ci * can.lai / DIM_LAYER)) / (OPTICS.ko * can.ci * can.lai / DIM_LAYER);
    OPTICS.po .= exp.(can._x_bnds .* OPTICS.ko .* can.ci .* can.lai) .* _fac_o;
    OPTICS.ps .= exp.(can._x_bnds .* OPTICS.ks .* can.ci .* can.lai) .* _fac_s;
    OPTICS.p_sunlit .= (view(OPTICS.ps,1:DIM_LAYER) .+ view(OPTICS.ps,2:DIM_LAYER+1)) ./ 2;

    _dso = sqrt( tand(angles.sza) ^ 2 + tand(angles.vza) ^ 2 - 2 * tand(angles.sza) * tand(angles.vza) * cosd(angles.vaa - angles.saa) );
    @inline _pdf(x::FT) where {FT<:AbstractFloat} = (
        _Σk = OPTICS.ko + OPTICS.ks;
        _Πk = OPTICS.ko * OPTICS.ks;
        _cl = can.ci * can.lai;
        _α  = _dso / HOT_SPOT * 2 / _Σk;

        if _dso == 0
            return exp( (_Σk - sqrt(_Πk)) * _cl * x )
        end;

        return exp( _Σk * _cl * x + sqrt(_Πk) * _cl / _α * (1 - exp(_α * x)) )
    );

    for _i in eachindex(can._x_bnds)
        OPTICS.pso[_i] = quadgk(_pdf, can._x_bnds[_i] - FT(1)/DIM_LAYER, can._x_bnds[_i]; rtol = 1e-2)[1] * DIM_LAYER;
    end;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-14: add method to use broadband PAR and NIR soil albedo for canopy_optical_properties!
#     2022-Jun-14: add method to use hyperspectral soil albedo for canopy_optical_properties!
#     2023-Mar-16: update both dd and sd from broadband values
#
#######################################################################################################################################################################################################
"""

    canopy_optical_properties!(can::HyperspectralMLCanopy{FT}, albedo::BroadbandSoilAlbedo{FT}) where {FT<:AbstractFloat}
    canopy_optical_properties!(can::HyperspectralMLCanopy{FT}, albedo::HyperspectralSoilAlbedo{FT}) where {FT<:AbstractFloat}

Updates lower soil boundary reflectance, given
- `can` `HyperspectralMLCanopy` type struct
- `albedo` `BroadbandSoilAlbedo` or `HyperspectralSoilAlbedo` type soil albedo

"""
canopy_optical_properties!(can::HyperspectralMLCanopy{FT}, albedo::BroadbandSoilAlbedo{FT}) where {FT<:AbstractFloat} = (
    (; OPTICS, WLSET) = can;

    OPTICS.ρ_dd[WLSET.IΛ_PAR,end] .= albedo.ρ_sw[1];
    OPTICS.ρ_dd[WLSET.IΛ_NIR,end] .= albedo.ρ_sw[2];
    OPTICS.ρ_sd[WLSET.IΛ_PAR,end] .= albedo.ρ_sw[1];
    OPTICS.ρ_sd[WLSET.IΛ_NIR,end] .= albedo.ρ_sw[2];

    return nothing
);

canopy_optical_properties!(can::HyperspectralMLCanopy{FT}, albedo::HyperspectralSoilAlbedo{FT}) where {FT<:AbstractFloat} = (
    (; OPTICS) = can;

    OPTICS.ρ_dd[:,end] .= albedo.ρ_sw;
    OPTICS.ρ_sd[:,end] .= albedo.ρ_sw;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-08: migrate the function from CanopyLayers
#     2022-Jun-09: rename the function from canopy_matrices! to canopy_optical_properties!
#     2022-Jun-09: move part of the short_wave! code into canopy_optical_properties!
#     2022-Jun-10: add function to compute longwave reflectance, transmittance, and emissivity
#     2022-Jun-29: use Leaves2D for the hyperspectral RT
#     2023-Mar-11: add code to account for the case of LAI == 0
#
#######################################################################################################################################################################################################
"""

    canopy_optical_properties!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaves2D{FT}}, soil::Soil{FT}) where {FT<:AbstractFloat}

Updates canopy optical properties (scattering coefficient matrices), given
- `can` `HyperspectralMLCanopy` type struct
- `leaves` Vector of `Leaves2D`
- `soil` Bottom soil boundary layer

"""
canopy_optical_properties!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaves2D{FT}}, soil::Soil{FT}) where {FT<:AbstractFloat} = (
    (; DIM_LAYER, OPTICS) = can;
    (; ALBEDO) = soil;
    @assert length(leaves) == DIM_LAYER "Number of leaves must be equal to the canopy layers!";

    if can.lai == 0
        OPTICS.ρ_dd .= 0;
        OPTICS.ρ_lw .= 0;
        OPTICS.ρ_sd .= 0;
        OPTICS.τ_dd .= 0;
        OPTICS.τ_lw .= 0;
        OPTICS.τ_sd .= 0;
        OPTICS._τ_ss = 0;
        canopy_optical_properties!(can, ALBEDO);
        OPTICS.ρ_lw[end] = ALBEDO.ρ_LW;

        return nothing
    end;

    _ilai = can.lai * can.ci / DIM_LAYER;

    # 1. update the scattering coefficients for different layers
    for _i in eachindex(leaves)
        OPTICS.σ_ddb[:,_i] .= OPTICS.ddb * leaves[_i].BIO.ρ_sw .+ OPTICS.ddf * leaves[_i].BIO.τ_sw;
        OPTICS.σ_ddf[:,_i] .= OPTICS.ddf * leaves[_i].BIO.ρ_sw .+ OPTICS.ddb * leaves[_i].BIO.τ_sw;
        OPTICS.σ_sdb[:,_i] .= OPTICS.sdb * leaves[_i].BIO.ρ_sw .+ OPTICS.sdf * leaves[_i].BIO.τ_sw;
        OPTICS.σ_sdf[:,_i] .= OPTICS.sdf * leaves[_i].BIO.ρ_sw .+ OPTICS.sdb * leaves[_i].BIO.τ_sw;
        OPTICS.σ_dob[:,_i] .= OPTICS.dob * leaves[_i].BIO.ρ_sw .+ OPTICS.dof * leaves[_i].BIO.τ_sw;
        OPTICS.σ_dof[:,_i] .= OPTICS.dof * leaves[_i].BIO.ρ_sw .+ OPTICS.dob * leaves[_i].BIO.τ_sw;
        OPTICS.σ_so[:,_i]  .= OPTICS.sob * leaves[_i].BIO.ρ_sw .+ OPTICS.sof * leaves[_i].BIO.τ_sw;
    end;

    # 2. update the transmittance and reflectance for single directions per layer (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    OPTICS._τ_ss  = exp(-1 * OPTICS.ks * _ilai);
    OPTICS._τ_dd .= exp.(-1 .* (1 .- OPTICS.σ_ddf) .* _ilai);
    OPTICS._τ_sd .= OPTICS.σ_sdf .* _ilai;  # TODO: use exp for these!
    OPTICS._ρ_dd .= OPTICS.σ_ddb .* _ilai;  # TODO: use exp for these!
    OPTICS._ρ_sd .= OPTICS.σ_sdb .* _ilai;  # TODO: use exp for these!

    # 3. update the effective reflectance per layer
    canopy_optical_properties!(can, ALBEDO);

    for _i in DIM_LAYER:-1:1
        _r_dd__ = view(OPTICS._ρ_dd,:,_i  );    # reflectance without correction
        _r_dd_i = view(OPTICS.ρ_dd ,:,_i  );    # reflectance of the upper boundary (i)
        _r_dd_j = view(OPTICS.ρ_dd ,:,_i+1);    # reflectance of the lower boundary (i+1)
        _r_sd__ = view(OPTICS._ρ_sd,:,_i  );    # reflectance without correction
        _r_sd_i = view(OPTICS.ρ_sd ,:,_i  );    # reflectance of the upper boundary (i)
        _r_sd_j = view(OPTICS.ρ_sd ,:,_i+1);    # reflectance of the lower boundary (i+1)
        _t_dd__ = view(OPTICS._τ_dd,:,_i  );    # transmittance without correction
        _t_dd_i = view(OPTICS.τ_dd ,:,_i  );    # transmittance of the layer (i)
        _t_sd__ = view(OPTICS._τ_sd,:,_i  );    # transmittance without correction
        _t_sd_i = view(OPTICS.τ_sd ,:,_i  );    # transmittance of the layer (i)
        _t_ss__ = OPTICS._τ_ss;                 # transmittance for directional->directional

        OPTICS._tmp_vec_λ .= 1 .- _r_dd__ .* _r_dd_j;

        _t_sd_i .= (_t_sd__ .+ _t_ss__ .* _r_sd_j .* _r_dd__) ./ OPTICS._tmp_vec_λ;             # it + it-jr-ir, then rescale
        _t_dd_i .= _t_dd__ ./ OPTICS._tmp_vec_λ;                                                # it, then rescale
        _r_sd_i .= _r_sd__ .+ _t_ss__ .* _r_sd_j .* _t_dd__ .+ _t_sd_i .* _r_dd_j .* _t_dd__;   # ir + it-jr-it(v) + it-jr_dd-it
        _r_dd_i .= _r_dd__ .+ _t_dd__ .* _r_dd_j .* _t_dd_i;                                    # ir + it-jr-it
    end;

    # 4. compute longwave effective emissivity, reflectance, and transmittance (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    for _i in eachindex(leaves)
        _σ_lw_b = OPTICS.ddb * leaves[_i].BIO.ρ_LW + OPTICS.ddf * leaves[_i].BIO.τ_LW;
        _σ_lw_f = OPTICS.ddf * leaves[_i].BIO.ρ_LW + OPTICS.ddb * leaves[_i].BIO.τ_LW;
        OPTICS._τ_lw[_i] = exp(-1 * (1 - _σ_lw_f) * _ilai);
        OPTICS._ρ_lw[_i] = _σ_lw_b * _ilai;
        OPTICS.ϵ[_i] = 1 - OPTICS._τ_lw[_i] - OPTICS._ρ_lw[_i];
    end;

    # 5. update the effective longwave reflectance and transmittance
    OPTICS.ρ_lw[end] = ALBEDO.ρ_LW;

    for _i in DIM_LAYER:-1:1
        _dnorm = 1 - OPTICS._ρ_lw[_i] * OPTICS.ρ_lw[_i+1];

        OPTICS.τ_lw[_i] = OPTICS._τ_lw[_i] / _dnorm;                                                    # it, then rescale
        OPTICS.ρ_lw[_i] = OPTICS._ρ_lw[_i] + OPTICS._τ_lw[_i] * OPTICS.ρ_lw[_i+1] * OPTICS.τ_lw[_i];    # ir + it-jr-it
    end;

    return nothing
);
