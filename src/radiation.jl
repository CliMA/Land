#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-22: add function to compute leaf level PAR and APAR
#     2022-Jan-13: use LeafBiophysics directly in the function rather than Leaf
#     2022-Feb-02: fix documentation
#     2022-Feb-02: unpack CONSTANTS only
#     2022-Jun-15: rename LeafBiophysics to HyperspectralLeafBiophysics to be more descriptive
#     2022-Jun-27: refactor the function to return PAR, APAR, and PPAR
#
#######################################################################################################################################################################################################
"""

    leaf_PAR(bio::HyperspectralLeafBiophysics{FT}, wls::WaveLengthSet{FT}, rad::HyperspectralRadiation{FT}; APAR_car::Bool = true) where {FT<:AbstractFloat}

Return leaf level PAR, APAR, and PPAR, given
- `bio` `ClimaCache.HyperspectralLeafBiophysics` type struct that contains leaf biophysical parameters
- `wls` `ClimaCache.WaveLengthSet` type struct that contains wave length bins
- `rad` `ClimaCache.HyperspectralRadiation` type struct that contains incoming radiation information
- `APAR_car` If true (default), account carotenoid absorption as PPAR; otherwise, PPAR is only by chlorophyll

---
# Examples
```julia
wls = WaveLengthSet{Float64}();
bio = HyperspectralLeafBiophysics{Float64}(wls);
rad = HyperspectralRadiation{Float64}();
par,apar,ppar = leaf_PAR(bio, wls, rad);
par,apar,ppar = leaf_PAR(bio, wls, rad; APAR_car=false);
```
"""
function leaf_PAR(bio::HyperspectralLeafBiophysics{FT}, wls::WaveLengthSet{FT}, rad::HyperspectralRadiation{FT}; APAR_car::Bool = true) where {FT<:AbstractFloat}
    @unpack IΛ_PAR, ΔΛ_PAR, Λ_PAR = wls;

    # PPAR absorption feature (after APAR is computed)
    _α_ppar = (APAR_car ? view(bio.α_cabcar, IΛ_PAR) : view(bio.α_cab, IΛ_PAR));

    # PAR, APAR, and PPAR energy from direct and diffuse light
    _e_par_dir   = view(rad.e_direct , IΛ_PAR);
    _e_par_diff  = view(rad.e_diffuse, IΛ_PAR);
    _e_apar_dir  = view(bio.α_sw, IΛ_PAR) .* _e_par_dir;
    _e_apar_diff = view(bio.α_sw, IΛ_PAR) .* _e_par_diff;
    _e_ppar_dir  = _α_ppar .* _e_apar_dir;
    _e_ppar_diff = _α_ppar .* _e_apar_diff;

    # PAR, APAR, and PPAR photons from direct and diffuse light
    _par_dir   = photon.(Λ_PAR, _e_par_dir  );
    _par_diff  = photon.(Λ_PAR, _e_par_diff );
    _apar_dir  = photon.(Λ_PAR, _e_apar_dir );
    _apar_diff = photon.(Λ_PAR, _e_apar_diff);
    _ppar_dir  = photon.(Λ_PAR, _e_ppar_dir );
    _ppar_diff = photon.(Λ_PAR, _e_ppar_diff);

    # total PAR and APAR in μmol photons m⁻² s⁻¹
    _Σpar_dir   = numerical∫(_par_dir  , ΔΛ_PAR) * 1000;
    _Σpar_diff  = numerical∫(_par_diff , ΔΛ_PAR) * 1000;
    _Σapar_dir  = numerical∫(_apar_dir , ΔΛ_PAR) * 1000;
    _Σapar_diff = numerical∫(_apar_diff, ΔΛ_PAR) * 1000;
    _Σppar_dir  = numerical∫(_ppar_dir , ΔΛ_PAR) * 1000;
    _Σppar_diff = numerical∫(_ppar_diff, ΔΛ_PAR) * 1000;

    return _Σpar_dir + _Σpar_diff, _Σapar_dir + _Σapar_diff, _Σppar_dir + _Σppar_diff
end
