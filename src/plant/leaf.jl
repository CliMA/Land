#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jan-14: refactor the Leaf structure within BIO, PRC, PSM as fields
#     2022-Jan-24: add p_CO₂_s to the structure
#     2022-Jan-24: fix documentation
#     2022-Feb-07: moved FLM to PRC
#     2022-May-25: add new field HS
#     2022-May-25: add new field WIDTH
#     2022-Jun-13: use Union instead of Abstract... for type definition
# Bug fixes:
#     2022-Jan-24: add FT control to p_CO₂_i
# To do
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
    "[`AbstractLeafBiophysics`](@ref) type leaf biophysical parameters"
    BIO::HyperspectralLeafBiophysics{FT}
    "[`LeafHydraulics`](@ref) type leaf hydraulic system"
    HS::LeafHydraulics{FT}
    "[`AbstractReactionCenter`](@ref) type photosynthesis reaction center"
    PRC::Union{VJPReactionCenter{FT}, CytochromeReactionCenter{FT}}
    "[`AbstractPhotosynthesisModel`](@ref) type photosynthesis model"
    PSM::Union{C3VJPModel{FT}, C4VJPModel{FT}, C3CytochromeModel{FT}}
    "Leaf width"
    WIDTH::FT

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
#     2022-Feb-07: remove fluorescence model from Leaf struct
#     2022-Feb-11: set default APAR = 1000
#     2022-Feb-11: add colimit option in constructor to enable quick deployment of quadratic colimitation
#     2022-May-25: add leaf hydraulic system into the constructor
#     2022-May-31: add steady state mode option to input options
#     2022-May-25: add new field WIDTH
#
#######################################################################################################################################################################################################
"""

    Leaf{FT}(psm::String, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); colimit::Bool = false, ssm::Bool = true) where {FT<:AbstractFloat}

Constructor for `Leaf`, given
- `psm` Photosynthesis model type, must be `C3`, `C3Cytochrome`, or `C4`
- `wls` [`WaveLengthSet`](@ref) type structure that determines the dimensions of leaf parameters
- `ssm` Whether the flow rate is at steady state

---
# Examples
```julia
leaf_c3 = Leaf{Float64}("C3");
leaf_c4 = Leaf{Float64}("C4");
leaf_cy = Leaf{Float64}("C3Cytochrome");
leaf_c3 = Leaf{Float64}("C3"; colimit = true);
leaf_c4 = Leaf{Float64}("C4"; colimit = true);
leaf_cy = Leaf{Float64}("C3Cytochrome"; colimit = true);
wls = WaveLengthSet{FT}(collect(400:10:2500));
leaf_c3 = Leaf{Float64}("C3", wls);
leaf_c4 = Leaf{Float64}("C4", wls);
leaf_cy = Leaf{Float64}("C3Cytochrome", wls);
```
"""
Leaf{FT}(psm::String, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); colimit::Bool = false, ssm::Bool = true) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C3Cytochrome", "C4"] "Photosynthesis model ID must be C3, C4, or C3Cytochrome!";

    if psm == "C3"
        _prc = VJPReactionCenter{FT}();
        _psm = C3VJPModel{FT}(colimit = colimit);
    elseif psm == "C3Cytochrome"
        _prc = CytochromeReactionCenter{FT}();
        _psm = C3CytochromeModel{FT}(colimit = colimit);
    elseif psm == "C4"
        _prc = VJPReactionCenter{FT}();
        _psm = C4VJPModel{FT}(colimit = colimit);
    end;

    return Leaf{FT}(
                HyperspectralLeafBiophysics{FT}(wls),            # BIO
                LeafHydraulics{FT}(ssm = ssm),      # HS
                _prc,                               # PRC
                _psm,                               # PSM
                FT(0.05),                           # WIDTH
                1000,                               # apar
                0.01,                               # g_H₂O_s
                T_25(),                             # t
                0.01,                               # g_CO₂
                3.0,                                # g_CO₂_b
                20,                                 # p_CO₂_i
                40,                                 # p_CO₂_s
                saturation_vapor_pressure(T_25()),  # p_H₂O_sat
                0)                                  # _t
);
