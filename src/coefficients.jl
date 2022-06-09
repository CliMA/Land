#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-07: migrate the function from CanopyLayers
#     2022-Jun-07: rename the function to be more descriptive extinction_scattering_coefficients
#     2022-Jun-07: add relative azimuth angle control
#     2022-Jun-07: return _Co, _Cs, _So, _Ss as well
# Sources
#     Verhoef (1998) Theory of radiative transfer models applied in optical remote sensing of vegetation canopies. Chapter 7
#
#######################################################################################################################################################################################################
"""

    extinction_scattering_coefficients(sza::FT, vza::FT, raa::FT, lia::FT) where {FT<:AbstractFloat}

Return the extinction and scattering coefficients (extinction coefficients from solar and viewing directions, and scattering coefficients for backward and forward directions, and some sin and cos
    products: _Co, _Cs, _So, _Ss), given
- `sza` Solar zenith angle in `°`
- `vza` Viewing zenith angle in `°`
- `raa` Relative azimuth angle in `°`
- `lia` Leaf inclination angle in `°`
"""
function extinction_scattering_coefficients(sza::FT, vza::FT, raa::FT, lia::FT) where {FT<:AbstractFloat}
    _π = FT(pi);

    # 1. compute the extinction coefficients (_ks for direct solar flux Es, _ko for flux-equivalent radiance Eo)
    _Cs = cosd(lia) * cosd(sza);
    _Ss = sind(lia) * sind(sza);
    _Co = cosd(lia) * cosd(vza);
    _So = sind(lia) * sind(vza);
    _βs = (_Cs >= _Ss ? _π : acos(-_Cs/_Ss));
    _βo = (_Co >= _So ? _π : acos(-_Co/_So));
    _Ds = (_βs < pi ? _Ss : _Cs);
    _Do = (0 < _βo < pi ? _So : -_Co/cos(_βo));
    _ks = 2 / _π / cosd(sza) * (_Cs * (_βs - _π/2) + _Ss * sin(_βs));
    _ko = 2 / _π / cosd(vza) * (_Co * (_βo - _π/2) + _So * sin(_βo));

    # 2. compute the scattering coefficients
    # 2.1 compute the Δ and β angles
    _Δ₁ = abs(_βs - _βo);
    _Δ₂ = _π - abs(_βs + _βo - _π);

    _ψ = deg2rad( abs(raa - 360*round(raa/360)) );
    if _ψ <= _Δ₁
        _β₁,_β₂,_β₃ = _ψ,_Δ₁,_Δ₂;
    elseif _Δ₁ < _ψ < _Δ₂
        _β₁,_β₂,_β₃ = _Δ₁,_ψ,_Δ₂;
    else
        _β₁,_β₂,_β₃ = _Δ₁,_Δ₂,_ψ;
    end

    # 2.2 compute the scattering coefficients
    _so = cosd(sza) * cosd(vza);
    _T₁ = 2 * _Cs * _Co + _Ss * _So * cosd(raa);
    _T₂ = sin(_β₂) * (2 * _Ds * _Do + _Ss * _So * cos(_β₁) * cos(_β₃));
    _F₁ = ((_π - _β₂) * _T₁ + _T₂) / _so;
    _F₂ = (-_β₂ * _T₁ + _T₂) / _so;

    # 2.3 compute the area scattering coefficient fractions (_sb for backward and _sf for forward)
    _sb = _F₁ / (2 * _π);
    _sf = _F₂ / (2 * _π);

    return _ks, _ko, _sb, _sf, _Co, _Cs, _So, _Ss
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-07: add function to update the extinction and scattering coefficient cache within canopy
#     2022-Jun-07: update _Co, _Cs, _So, _Ss as well
#
#######################################################################################################################################################################################################
"""

    extinction_scattering_coefficients!(can::HyperspectralMLCanopy{FT}, angles::SunSensorGeometry{FT}) where {FT<:AbstractFloat}

Update the extinction and scattering coefficients, given
- `can` `HyperspectralMLCanopy` type canopy
- `angles` `SunSensorGeometry` type angles
"""
function extinction_scattering_coefficients!(can::HyperspectralMLCanopy{FT}, angles::SunSensorGeometry{FT}) where {FT<:AbstractFloat}
    @unpack OPTICS, Θ_INCL = can;

    for _i in eachindex(Θ_INCL)
        OPTICS._ks[_i], OPTICS._ko[_i], OPTICS._sb[_i], OPTICS._sf[_i], OPTICS._Co[_i], OPTICS._Cs[_i], OPTICS._So[_i], OPTICS._Ss[_i] =
            extinction_scattering_coefficients(angles.sza, angles.vza, angles.vaa - angles.saa, Θ_INCL[_i]);
    end;

    return nothing
end
