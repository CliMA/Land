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
# Changes to this struct
# General
#     2022-Jun-15: add struct for broadband leaf biophysics
#     2022-Jun-24: add leaf emissivity constant
#     2022-Jul-19: add field lma
#     2022-Jul-20: use α and ϵ instead of upper case greek Α and Ε
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf biophysical traits used to run leaf reflectance and transmittance.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BroadbandLeafBiophysics{FT<:AbstractFloat} <: AbstractLeafBiophysics{FT}
    # General information of leaf biophysics
    "Broadband absorption fraction at the NIR region"
    α_NIR::FT = 0.2
    "Broadband absorption fraction at the PAR region"
    α_PAR::FT = 0.8
    "Emissivity for longwave radiation"
    ϵ_LW::FT = 0.97

    # Prognostic variables
    "Dry matter content (dry leaf mass per unit area) `[g cm⁻²]`"
    lma::FT = 0.012
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Aug-04: refactor the structure with constants, variables, and temporary cache
#     2021-Aug-04: add concentrations and characteristic curves altogether
#     2021-Aug-10: add CBC and PRO supoort
#     2021-Sep-30: rename LeafBio to LeafBiophysics to be more specific
#     2021-Oct-21: rename f_sense and K_SENES to brown and K_BROWN
#     2021-Nov-24: tease apart the characteristic absorption curves to HyperspectralAbsorption
#     2022-Jun-15: rename struct to HyperspectralLeafBiophysics to distinguish from BroadbandLeafBiophysics
#     2022-Jul-19: add dimension control to struct
#     2022-Jul-20: remove field: l_H₂O (unit cm)
#     2022-Jul-20: rename ρ_lw and τ_lw to ρ_LW and τ_LW
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf biophysical traits used to run leaf reflectance and transmittance.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct HyperspectralLeafBiophysics{FT<:AbstractFloat} <: AbstractLeafBiophysics{FT}
    # Dimensions
    "Dimension of SIF wave length bins"
    DIM_SIF::Int = 29
    "Dimension of SIF excitation wave length bins"
    DIM_SIFE::Int = 45
    "Dimension of short wave length bins"
    DIM_WL::Int = 114

    # General information of leaf biophysics
    "Leaf mesophyll structural parameter that describes mesophyll reflectance and transmittance"
    MESOPHYLL_N::FT = 1.4
    "Doubling adding layers"
    NDUB::Int = 10
    "Broadband thermal reflectance, related to blackbody emittance `[-]`"
    ρ_LW::FT = 0.01
    "Broadband thermal transmission, related to blackbody emittance `[-]`"
    τ_LW::FT = 0.01

    # Prognostic variables
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
    "Dry matter content (dry leaf mass per unit area) `[g cm⁻²]`"
    lma::FT = 0.012
    "Protein content in lma (pro = lma - cbc) `[g cm⁻²]`"
    pro::FT = 0

    # Diagnostic variables
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
    "Shortwave leaf reflectance `[-]`"
    ρ_sw::Vector{FT} = zeros(FT, DIM_WL)
    "Shortwave leaf transmission `[-]`"
    τ_sw::Vector{FT} = zeros(FT, DIM_WL)
end
