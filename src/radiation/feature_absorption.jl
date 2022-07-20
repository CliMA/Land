#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Aug-04: refactor the structure with constants, variables, and temporary cache
#     2021-Aug-04: add concentrations and characteristic curves altogether
#     2021-Aug-10: add CBC and PRO supoort
#     2021-Agu-10: add constructors within the structure rather than initialize it externally
#     2021-Sep-30: rename LeafBio to LeafBiophysics to be more specific
#     2021-Oct-19: sort variable to prognostic and dignostic catergories
#     2021-Oct-21: rename f_sense and K_SENES to brown and K_BROWN
#     2021-Nov-24: tease apart the characteristic absorption curves to HyperspectralAbsorption
#     2022-Jul-20: use kwdef for the constructor
#     2022-Jul-20: add field DATASET to struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Immutable struct that contains leaf biophysical traits used to run leaf reflection and transmittance.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct HyperspectralAbsorption{FT<:AbstractFloat}
    # File path to the Netcdf dataset
    "File path to the Netcdf dataset"
    DATASET::String = LAND_2021

    # parameters that do not change with time
    "Specific absorption coefficients of anthocynanin `[-]`"
    K_ANT::Vector{FT} = read_nc(DATASET, "K_ANT")
    "Specific absorption coefficients of senescent material (brown pigments) `[-]`"
    K_BROWN::Vector{FT} = read_nc(DATASET, "K_BROWN")
    "Specific absorption coefficients of chlorophyll a and b `[-]`"
    K_CAB::Vector{FT} = read_nc(DATASET, "K_CAB")
    "Specific absorption coefficients of violaxanthin carotenoid `[-]`"
    K_CAR_V::Vector{FT} = read_nc(DATASET, "K_CAR_V")
    "Specific absorption coefficients of zeaxanthin carotenoid `[-]`"
    K_CAR_Z::Vector{FT} = read_nc(DATASET, "K_CAR_Z")
    "Specific absorption coefficients of carbon-based constituents `[-]`"
    K_CBC::Vector{FT} = read_nc(DATASET, "K_CBC")
    "Specific absorption coefficients of water `[-]`"
    K_H₂O::Vector{FT} = read_nc(DATASET, "K_H₂O")
    "Specific absorption coefficients of dry matter `[-]`"
    K_LMA::Vector{FT} = read_nc(DATASET, "K_LMA")
    "Specific absorption coefficients of protein `[-]`"
    K_PRO::Vector{FT} = read_nc(DATASET, "K_PRO")
    "Specific absorption coefficients of PS I and II `[-]`"
    K_PS::Vector{FT} = read_nc(DATASET, "K_PS")
    "Refractive index `[-]`"
    NR::Vector{FT} = read_nc(DATASET, "NR")
end
