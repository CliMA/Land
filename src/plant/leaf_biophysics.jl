#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-15: add abstract type for leaf biophysics
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractSoilAlbedo:
- [`BroadbandLeafBiophysics`](@ref)
- [`HyperspectralLeafBiophysics`](@ref)

"""
abstract type AbstractLeafBiophysics{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-15: add struct for broadband leaf biophysics
#     2022-Jun-24: add leaf emissivity constant
#     2022-Jul-18: use kwdef for the constructor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf biophysical traits used to run leaf reflectance and transmittance.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BroadbandLeafBiophysics{FT} <: AbstractLeafBiophysics{FT}
    # parameters that do not change with time
    "Broadband absorption fraction at the NIR region"
    Α_NIR::FT = 0.2
    "Broadband absorption fraction at the PAR region"
    Α_PAR::FT = 0.8
    "Emissivity for longwave radiation"
    Ε_LW::FT = 0.97
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2021-Aug-04: refactor the structure with constants, variables, and temporary cache
#     2021-Aug-04: add concentrations and characteristic curves altogether
#     2021-Aug-10: add CBC and PRO supoort
#     2021-Agu-10: add constructors within the structure rather than initialize it externally
#     2021-Sep-30: rename LeafBio to LeafBiophysics to be more specific
#     2021-Oct-19: sort variable to prognostic and dignostic catergories
#     2021-Oct-21: rename f_sense and K_SENES to brown and K_BROWN
#     2021-Nov-24: tease apart the characteristic absorption curves to HyperspectralAbsorption
#     2022-Jun-15: rename struct to HyperspectralLeafBiophysics to distinguish from BroadbandLeafBiophysics
#     2022-Jul-18: use kwdef for the constructor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf biophysical traits used to run leaf reflectance and transmittance.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct HyperspectralLeafBiophysics{FT} <: AbstractLeafBiophysics{FT}
    # parameters that do not change with time
    "Leaf mesophyll structural parameter that describes mesophyll reflectance and transmittance"
    MESOPHYLL_N::FT = 1.4
    "Doubling adding layers"
    NDUB::Int = 10

    # prognostic variables that change with time
    "Anthocyanin content `[μg cm⁻²]`"
    ant::FT = 0
    "Senescent material (brown pigments) fraction `[-]`"
    brown::FT = 0
    "Chlorophyll a and b content `[μg cm⁻²]`"
    cab::FT = 40
    "Carotenoid content `[μg cm⁻²]`"
    car::FT = 40 / 7
    "Carbon-based constituents in lma `[g cm⁻²]`"
    cbc::FT = 0
    "Zeaxanthin fraction in Carotenoid (1=all Zeaxanthin, 0=all Violaxanthin) `[-]`"
    f_zeax::FT = 0
    "Equivalent water thickness `[cm]`" # TODO: remove this and use LeafHydraulics field for the calculation
    l_H₂O::FT = 0.01
    "Dry matter content (dry leaf mass per unit area) `[g cm⁻²]`"
    lma::FT = 0.012
    "Protein content in lma (pro = lma - cbc) `[g cm⁻²]`"
    pro::FT = 0

    # dignostic variables that change with time
    "Specific absorption coefficients of all materials"
    k_all::Vector{FT}
    "Fluorescence excitation matrix backwards `[-]`"
    mat_b::Matrix{FT}
    "Fluorescence excitation matrix forwards `[-]`"
    mat_f::Matrix{FT}
    "Relative absorption by Chlorophyll `[-]`"
    α_cab::Vector{FT}
    "Relative absorption by Chlorophyll+Carotenoid `[-]`"
    α_cabcar::Vector{FT}
    "Shortwave absorption, 1 .- ρ_sw .- τ_sw  `[-]`"
    α_sw::Vector{FT}
    "Broadband thermal reflectance, related to blackbody emittance `[-]`"
    ρ_lw::FT = 0.01
    "Shortwave leaf reflectance `[-]`"
    ρ_sw::Vector{FT}
    "Broadband thermal transmission, related to blackbody emittance `[-]`"
    τ_lw::FT = 0.01
    "Shortwave leaf transmission `[-]`"
    τ_sw::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2021-Nov-24: migrate the constructor from CanopyLayers
#     2022-Jun-15: rename struct to HyperspectralLeafBiophysics to distinguish from BroadbandLeafBiophysics
#     2022-Jul-18: use external constructor for those 2D matrices based on WaveLengthSet
#
#######################################################################################################################################################################################################
"""

    HyperspectralLeafBiophysics{FT}(wls::WaveLengthSet{FT} = WaveLengthSet{FT}()) where {FT<:AbstractFloat}

Constructor for `HyperspectralLeafBiophysics`, given
- `wls` [`WaveLengthSet`](@ref) type structure

"""
HyperspectralLeafBiophysics{FT}(wls::WaveLengthSet{FT} = WaveLengthSet{FT}()) where {FT<:AbstractFloat} = (
    @unpack NΛ, NΛ_SIF, NΛ_SIFE = wls;

    return HyperspectralLeafBiophysics{FT}(
                k_all    = zeros(FT, NΛ),
                mat_b    = zeros(FT, NΛ_SIF, NΛ_SIFE),
                mat_f    = zeros(FT, NΛ_SIF, NΛ_SIFE),
                α_cab    = zeros(FT, NΛ),
                α_cabcar = zeros(FT, NΛ),
                α_sw     = zeros(FT, NΛ),
                ρ_sw     = zeros(FT, NΛ),
                τ_sw     = zeros(FT, NΛ)
    )
);
