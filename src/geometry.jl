#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-07: migrate the function from CanopyLayers
#     2022-Jun-07: clean the function
#     2022-Jun-08: add documentation
#
#######################################################################################################################################################################################################
"""

    canopy_geometry!(can::HyperspectralMLCanopy{FT}, angles::SunSensorGeometry{FT}) where {FT<:AbstractFloat}

Updates canopy optical properties (extinction coefficients for direct and diffuse light) based on the SAIL model, given
- `can` `HyperspectralMLCanopy` type struct
- `angles` `SunSensorGeometry` type struct
"""
function canopy_geometry!(can::HyperspectralMLCanopy{FT}, angles::SunSensorGeometry{FT}) where {FT<:AbstractFloat}
    @unpack HOT_SPOT, N_LAYER, OPTICS, P_INCL, Θ_AZI = can;

    # 1. update the canopy optical properties related to extinction and scattering coefficients
    extinction_scattering_coefficients!(can, angles);

    OPTICS.ko  = P_INCL' * OPTICS._ko;
    OPTICS.ks  = P_INCL' * OPTICS._ks;
    OPTICS.sob = P_INCL' * OPTICS._sb
    OPTICS.sof = P_INCL' * OPTICS._sf
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

    OPTICS._fs_fo .= OPTICS.fs .* OPTICS.fo;
    OPTICS._abs_fs_fo .= abs.(OPTICS._fs_fo);

    # 3. update the viewing fraction ps, po, and pso
    _fac_s = (1 - exp(-OPTICS.ks * can.ci * can.lai / N_LAYER)) / (OPTICS.ks * can.ci * can.lai / N_LAYER);
    _fac_o = (1 - exp(-OPTICS.ko * can.ci * can.lai / N_LAYER)) / (OPTICS.ko * can.ci * can.lai / N_LAYER);
    OPTICS.po .= exp.(can._x_bnds * OPTICS.ko * can.ci * can.lai) * _fac_o;
    OPTICS.ps .= exp.(can._x_bnds * OPTICS.ks * can.ci * can.lai) * _fac_s;

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
        OPTICS.pso[_i] = quadgk(_pdf, can._x_bnds[_i] - FT(1)/N_LAYER, can._x_bnds[_i]; rtol = 1e-2)[1] * N_LAYER;
    end;

    return nothing
end
