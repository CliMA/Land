#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: rename the function from leaf_ETR! to to photosystem_electron_transport! to be more specific
#
#######################################################################################################################################################################################################
"""
This function updates the electron transport rates in photosynthesis reaction center, and thus to calculate photosynthetic rate. Supported methods are

$(METHODLIST)

"""
function photosystem_electron_transport! end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: use apar as an input rather than field from Leaf, and this allows for more modular operations
#     2022-Jan-14: remove examples from the doc, as this function is not meant to be public
#     2022-Jan-14: unpack CONSTANT from the input variables only
#
#######################################################################################################################################################################################################
"""
    photosystem_electron_transport!(psm::C3VJPModel{FT}, rc::PhotosynthesisReactionCenter{FT}, apar::FT) where {FT<:AbstractFloat}

Update the electron transport rates, given
- `psm` `C3VJPModel` type C3 photosynthesis model
- `rc` Photosynthesis system reaction center
- `apar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
"""
photosystem_electron_transport!(psm::C3VJPModel{FT}, rc::PhotosynthesisReactionCenter{FT}, apar::FT) where {FT<:AbstractFloat} = (
    @unpack F_PSII, Φ_PSII_MAX = rc;

    psm.j_pot = F_PSII * Φ_PSII_MAX * apar;
    psm.j = min(psm.j_pot, psm.j_max);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: use apar as an input rather than field from Leaf, and this allows for more modular operations
#     2022-Jan-14: remove examples from the doc, as this function is not meant to be public
#     2022-Jan-14: unpack CONSTANT from the input variables only
#
#######################################################################################################################################################################################################
"""
    photosystem_electron_transport!(psm::C4VJPModel{FT}, rc::PhotosynthesisReactionCenter{FT}, apar::FT) where {FT<:AbstractFloat}

Update the electron transport rates, given
- `psm` `C4VJPModel` type C4 photosynthesis model
- `rc` Photosynthesis system reaction center
- `apar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
"""
photosystem_electron_transport!(psm::C4VJPModel{FT}, rc::PhotosynthesisReactionCenter{FT}, apar::FT) where {FT<:AbstractFloat} = (
    @unpack F_PSII, Φ_PSII_MAX = rc;

    psm.j_pot = F_PSII * Φ_PSII_MAX * apar;
    psm.j = psm.j_pot;

    return nothing
);
