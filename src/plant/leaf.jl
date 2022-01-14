#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2021-Jan-14: refactor the Leaf structure within BIO, PRC, PSM as fields
# To do
#     TODO: add leaf physiological parameters as a field well
#     TODO: add leaf hydraulics as a field as well
#     TODO: link leaf water content to BIO_PHYSICS.l_H₂O
#
#######################################################################################################################################################################################################
"""
$(TYPEDEF)

Structure to save leaf parameters

# Fields
$(TYPEDFIELDS)
"""
mutable struct Leaf{FT<:AbstractFloat}
    # parameters that do not change with time
    "leaf biophysical parameters"
    BIO::LeafBiophysics{FT}
    "Photosynthesis reaction center"
    PRC::PhotosynthesisReactionCenter{FT}
    "Photosynthesis model"
    PSM::AbstractPhotosynthesisModel{FT}

    # prognostic variables that change with time
    "Current leaf temperature"
    t::FT

    # dignostic variables that change with time
    "Total leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂::FT
    "Leaf internal CO₂ partial pressure `[Pa]`"
    p_CO₂_i
    "Saturation H₂O vapor pressure, need to update with temperature and leaf water pressure `[Pa]`"
    p_H₂O_sat::FT

    # caches to speed up calculations
    "Last leaf temperature. If different from t, then make temperature correction"
    _t::FT
end


Leaf{FT}(psm::String, wls::WaveLengthSet{FT} = WaveLengthSet{FT}()) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C4"];

    _bio = LeafBiophysics{FT}(wls);
    _prc = PhotosynthesisReactionCenter{FT}();
    _t   = T_25();
    _p   = saturation_vapor_pressure(T_25());

    if psm == "C3"
        return Leaf{FT}(_bio, _prc, C3VJPModel{FT}(), _t, 0, 0, _p, 0)
    end;

    if psm == "C4"
        return Leaf{FT}(_bio, _prc, C4VJPModel{FT}(), _t, 0, 0, _p, 0)
    end;
);
