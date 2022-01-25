#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jan-14: refactor the Leaf structure within BIO, PRC, PSM as fields
#     2022-Jan-24: add p_CO₂_s to the structure
#     2022-Jan-24: fix documentation
# Bug fixes:
#     2022-Jan-24: add FT control to p_CO₂_i
# To do
#     TODO: add leaf physiological parameters as a field well
#     TODO: add leaf hydraulics as a field as well
#     TODO: link leaf water content to BIO_PHYSICS.l_H₂O
#     TODO: add auto reference in fields
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
    "[`LeafBiophysics`](@ref) type leaf biophysical parameters"
    BIO::LeafBiophysics{FT}
    "[`AbstractFluorescenceModel`](@ref) type leaf fluorescence model"
    FLM::AbstractFluorescenceModel{FT}
    "Photosynthesis reaction center"
    PRC::AbstractReactionCenter{FT}
    "Photosynthesis model"
    PSM::AbstractPhotosynthesisModel{FT}

    # prognostic variables that change with time
    "Absorbed photosynthetically active radiation `[μmol m⁻² s⁻¹]`"
    apar::FT
    "Stomatal conductance to water vapor `[mol m⁻² s⁻¹]`"
    g_H₂O_s::FT
    "Current leaf temperature"
    t::FT

    # dignostic variables that change with time
    "Total leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂::FT
    "Boundary leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂_b::FT
    "Leaf internal CO₂ partial pressure `[Pa]`"
    p_CO₂_i::FT
    "Leaf surface CO₂ partial pressure `[Pa]`"
    p_CO₂_s::FT
    "Saturation H₂O vapor pressure, need to update with temperature and leaf water pressure `[Pa]`"
    p_H₂O_sat::FT

    # caches to speed up calculations
    "Last leaf temperature. If different from t, then make temperature correction"
    _t::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jan-14: add C3 and C4 constructors
#     2022-Jan-24: add C3Cytochrome constructor
#     2022-Jan-24: add p_CO₂_s to the constructor
#     2022-Jan-24: add documentation
#
#######################################################################################################################################################################################################
"""

    Leaf{FT}(psm::String, wls::WaveLengthSet{FT} = WaveLengthSet{FT}()) where {FT<:AbstractFloat}

Constructor for `Leaf`, given
- `psm` Photosynthesis model type, must be `C3`, `C3Cytochrome`, or `C4`
- `wls` [`WaveLengthSet`](@ref) type structure that determines the dimensions of leaf parameters

---
# Examples
```julia
leaf_c3 = Leaf{Float64}("C3");
leaf_c4 = Leaf{Float64}("C4");
leaf_cy = Leaf{Float64}("C3Cytochrome");
wls = WaveLengthSet{FT}(collect(400:10:2500));
leaf_c3 = Leaf{Float64}("C3", wls);
leaf_c4 = Leaf{Float64}("C4", wls);
leaf_cy = Leaf{Float64}("C3Cytochrome", wls);
```
"""
Leaf{FT}(psm::String, wls::WaveLengthSet{FT} = WaveLengthSet{FT}()) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C3Cytochrome", "C4"];

    _bio = LeafBiophysics{FT}(wls);
    _t   = T_25();
    _p   = saturation_vapor_pressure(_t);

    _g_lc = 0.01;
    _g_bc = 3.0;
    _p_i  = 20.0;
    _p_s  = 40.0;

    if psm == "C3"
        return Leaf{FT}(_bio, VanDerTolFluorescenceModel{FT}(), VJPReactionCenter{FT}(), C3VJPModel{FT}(), 0, 0, _t, _g_lc, _g_bc, _p_i, _p_s, _p, 0)
    end;

    if psm == "C3Cytochrome"
        return Leaf{FT}(_bio, CytochromeFluorescenceModel{FT}(), CytochromeReactionCenter{FT}(), C3CytochromeModel{FT}(), 0, 0, _t, _g_lc, _g_bc, _p_i, _p_s, _p, 0)
    end;

    if psm == "C4"
        return Leaf{FT}(_bio, VanDerTolFluorescenceModel{FT}(), VJPReactionCenter{FT}(), C4VJPModel{FT}(), 0, 0, _t, _g_lc, _g_bc, _p_i, _p_s, _p, 0)
    end;
);
