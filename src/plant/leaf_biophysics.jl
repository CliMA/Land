#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-15: add abstract type for leaf biophysics
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractLeafBiophysics:
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
#     2022-Jul-19: add field lma
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf biophysical traits used to run leaf reflectance and transmittance.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BroadbandLeafBiophysics{FT<:AbstractFloat} <: AbstractLeafBiophysics{FT}
    # parameters that do not change with time
    "Broadband absorption fraction at the NIR region"
    Α_NIR::FT = 0.2
    "Broadband absorption fraction at the PAR region"
    Α_PAR::FT = 0.8
    "Emissivity for longwave radiation"
    Ε_LW::FT = 0.97

    # prognostic variables that change with time
    "Dry matter content (dry leaf mass per unit area) `[g cm⁻²]`"
    lma::FT = 0.012
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
#     2022-Jul-19: add dimension control to struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf biophysical traits used to run leaf reflectance and transmittance.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct HyperspectralLeafBiophysics{FT<:AbstractFloat} <: AbstractLeafBiophysics{FT}
    # dimensions
    "Dimension of SIF wave length bins"
    DIM_SIF::Int = 29
    "Dimension of SIF excitation wave length bins"
    DIM_SIFE::Int = 45
    "Dimension of short wave length bins"
    DIM_WL::Int = 114

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
    k_all::Vector{FT} = zeros(FT, DIM_WL)
    "Fluorescence excitation matrix backwards `[-]`"
    mat_b::Matrix{FT} = zeros(FT, DIM_SIF, DIM_SIFE)
    "Fluorescence excitation matrix forwards `[-]`"
    mat_f::Matrix{FT} = zeros(FT, DIM_SIF, DIM_SIFE)
    "Relative absorption by Chlorophyll `[-]`"
    α_cab::Vector{FT} = zeros(FT, DIM_WL)
    "Relative absorption by Chlorophyll+Carotenoid `[-]`"
    α_cabcar::Vector{FT} = zeros(FT, DIM_WL)
    "Shortwave absorption, 1 .- ρ_sw .- τ_sw  `[-]`"
    α_sw::Vector{FT} = zeros(FT, DIM_WL)
    "Broadband thermal reflectance, related to blackbody emittance `[-]`"
    ρ_lw::FT = 0.01
    "Shortwave leaf reflectance `[-]`"
    ρ_sw::Vector{FT} = zeros(FT, DIM_WL)
    "Broadband thermal transmission, related to blackbody emittance `[-]`"
    τ_lw::FT = 0.01
    "Shortwave leaf transmission `[-]`"
    τ_sw::Vector{FT} = zeros(FT, DIM_WL)
end
