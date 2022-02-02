#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-22: add function to compute leaf level PAR and APAR
#     2022-Jan-13: use LeafBiophysics directly in the function rather than Leaf
#     2022-Feb-02: fix documentation
#
#######################################################################################################################################################################################################
"""

    leaf_PAR(bio::LeafBiophysics{FT}, wls::WaveLengthSet{FT}, rad::HyperspectralRadiation{FT}; APAR_car::Bool = true) where {FT<:AbstractFloat}

Return leaf level PAR and APAR, given
- `bio` [`LeafBiophysics`](@ref) type struct that contains leaf biophysical parameters
- `wls` [`WaveLengthSet`](@ref) type struct that contains wave length bins
- `rad` [`HyperspectralRadiation`](@ref) type struct that contains incoming radiation information
- `APAR_car` If true (default), account carotenoid absorption as APAR; otherwise, APAR is only by chlorophyll

---
# Examples
```julia
wls = WaveLengthSet{Float64}();
bio = LeafBiophysics{Float64}(wls);
rad = HyperspectralRadiation{Float64}();
par,apar = leaf_PAR(bio, wls, rad);
par,apar = leaf_PAR(bio, wls, rad; APAR_car=false);
```
"""
function leaf_PAR(bio::LeafBiophysics{FT}, wls::WaveLengthSet{FT}, rad::HyperspectralRadiation{FT}; APAR_car::Bool = true) where {FT<:AbstractFloat}
    @unpack α_cab, α_cabcar, α_SW = bio;
    @unpack e_direct, e_diffuse = rad;
    @unpack IΛ_PAR, ΔΛ_PAR, Λ_PAR = wls;

    # APAR absorption feature
    _α = (APAR_car ? view(α_cabcar, IΛ_PAR) : view(α_cab, IΛ_PAR));

    # PAR energy from direct  and diffuse light
    _e_par_dir  = view(e_direct , IΛ_PAR) .* view(α_SW, IΛ_PAR);
    _e_par_diff = view(e_diffuse, IΛ_PAR) .* view(α_SW, IΛ_PAR);
    _par_dir  = photon.(Λ_PAR, _e_par_dir );
    _par_diff = photon.(Λ_PAR, _e_par_diff);

    # absorbed PAR energy from direct and diffuse light
    _apar_dir  = _α .* _par_dir;
    _apar_diff = _α .* _par_diff;

    # total PAR and APAR in μmol photons m⁻² s⁻¹
    _∑par_dir   = numerical∫(_par_dir, ΔΛ_PAR);
    _∑par_diff  = numerical∫(_par_diff, ΔΛ_PAR);
    _∑apar_dir  = numerical∫(_apar_dir, ΔΛ_PAR);
    _∑apar_diff = numerical∫(_apar_diff, ΔΛ_PAR);

    return 1000 * (_∑par_dir + _∑par_diff), 1000 * (_∑apar_dir + _∑apar_diff)
end
