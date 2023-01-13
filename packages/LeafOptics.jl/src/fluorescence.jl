#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Jul-08: add leaf level SIF simulation
#     2021-Jul-08: use mat_b and mat_f for SIF at backward and forward directions
#     2021-Aug-05: add option to sumulate SIF in photon to photon mode
#     2021-Oct-22: refactor the function to leaf_SIF to return the SIFs directly
#     2022-Jan-13: use LeafBiophysics directly in the function rather than Leaf
#     2022-Jun-15: rename LeafBiophysics to HyperspectralLeafBiophysics to be more descriptive
#
#######################################################################################################################################################################################################
"""

    leaf_SIF(bio::HyperspectralLeafBiophysics{FT}, wls::WaveLengthSet{FT}, rad::HyperspectralRadiation{FT}, ϕ::FT = FT(0.01); ϕ_photon::Bool = true) where {FT<:AbstractFloat}

Return the leaf level SIF at backward and forward directions, given
- `bio` `HyperspectralLeafBiophysics` type struct that contains leaf biophysical parameters
- `wls` `WaveLengthSet` type struct that contains wave length bins
- `rad` `HyperspectralRadiation` type struct that contains incoming radiation information
- `ϕ` Fluorescence quantum yield
- `ϕ_photon` If true (default), convert photon to photon when computing SIF; otherwise, convert energy to energy

---
# Examples
```julia
wls = EmeraldNamespace.WaveLengthSet{Float64}();
bio = EmeraldNamespace.HyperspectralLeafBiophysics{Float64}();
rad = EmeraldNamespace.HyperspectralRadiation{Float64}();
sif_b,sif_f = leaf_SIF(bio, wls, rad, 0.01);
sif_b,sif_f = leaf_SIF(bio, wls, rad, 0.01; ϕ_photon=false);
```

"""
function leaf_SIF(bio::HyperspectralLeafBiophysics{FT}, wls::WaveLengthSet{FT}, rad::HyperspectralRadiation{FT}, ϕ::FT = FT(0.01); ϕ_photon::Bool = true) where {FT<:AbstractFloat}
    (; IΛ_SIFE, ΔΛ_SIFE, Λ_SIF, Λ_SIFE) = wls;

    # calculate the excitation energy and photons
    _e_excitation = (view(rad.e_direct, IΛ_SIFE) .+ view(rad.e_diffuse, IΛ_SIFE)) .* ΔΛ_SIFE;

    # convert energy to energy using the matrices
    if !ϕ_photon
        _sif_b = bio.mat_b * _e_excitation * ϕ / FT(pi);
        _sif_f = bio.mat_f * _e_excitation * ϕ / FT(pi);

        return _sif_b, _sif_f
    end;

    # convert energy to photon
    _phot_excitation = photon.(Λ_SIFE, _e_excitation);

    # convert photon to photon using the matrices
    _phot_b = bio.mat_b * _phot_excitation * ϕ / FT(pi);
    _phot_f = bio.mat_f * _phot_excitation * ϕ / FT(pi);

    # convert photon to back to energy
    _sif_b = energy.(Λ_SIF, _phot_b);
    _sif_f = energy.(Λ_SIF, _phot_f);

    return _sif_b, _sif_f
end
