#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Aug-10: refactor the structure with renamed fields
#     2021-Aug-10: add a constructor within the structure to avoid external initialization
#     2021-Oct-19: sort variable to prognostic and dignostic catergories
#     2022-Jul-20: use kwdef for the constructor
#     2022-Jul-20: add field DATASET to struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Immutable structure that stores wave length information.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct WaveLengthSet{FT<:AbstractFloat}
    # File path to the Netcdf dataset
    "File path to the Netcdf dataset"
    DATASET::String = LAND_2021

    # Constants
    "Wavelength limits for NIR `[nm]`"
    WL_NIR::Vector{FT} = FT[700, 2500]
    "Wavelength limits for PAR `[nm]`"
    WL_PAR::Vector{FT} = FT[400, 750]
    "Wavelength limits for SIF emission `[nm]`"
    WL_SIF::Vector{FT} = FT[640, 850]
    "Wavelength limits for SIF excitation `[nm]`"
    WL_SIFE::Vector{FT} = FT[400, 750]
    "Wavelength (bins) `[nm]`"
    Λ::Vector{FT} = read_nc(DATASET, "WL")
    "Lower boundary wavelength `[nm]`"
    Λ_LOWER::Vector{FT} = read_nc(DATASET, "WL_LOWER")
    "Upper boundary wavelength `[nm]`"
    Λ_UPPER::Vector{FT} = read_nc(DATASET, "WL_UPPER")

    # Indices
    "Indicies of Λ_NIR in Λ"
    IΛ_NIR::Vector{Int} = findall( WL_NIR[1] .<= Λ .<= WL_NIR[2] )
    "Indicies of Λ_PAR in Λ"
    IΛ_PAR::Vector{Int} = findall( WL_PAR[1] .<= Λ .<= WL_PAR[2] )
    "Indicies of Λ_SIF in Λ"
    IΛ_SIF::Vector{Int} = findall( WL_SIF[1] .<= Λ .<= WL_SIF[2] )
    "Indicies of Λ_SIFE in Λ"
    IΛ_SIFE::Vector{Int} = findall( WL_SIFE[1] .<= Λ .<= WL_SIFE[2] )

    # Dimensions
    "Number of wavelength bins for NIR"
    DIM_NIR::Int = length(IΛ_NIR)
    "Number of wavelength bins for PAR"
    DIM_PAR::Int = length(IΛ_PAR)
    "Number of wavelength bins for SIF"
    DIM_SIF::Int = length(IΛ_SIF)
    "Number of wavelength bins for SIFE"
    DIM_SIFE::Int = length(IΛ_SIFE)
    "Number of wavelength bins"
    DIM_WL::Int = length(Λ)

    # Constants based on the ones above
    "Differential wavelength `[nm]`"
    ΔΛ::Vector{FT} = Λ_UPPER .- Λ_LOWER
    "Differential wavelength for PAR `[nm]`"
    ΔΛ_PAR::Vector{FT} = ΔΛ[IΛ_PAR]
    "Differential wavelength for SIF excitation `[nm]`"
    ΔΛ_SIFE::Vector{FT} = ΔΛ[IΛ_SIFE]
    "Wavelength bins for PAR `[nm]`"
    Λ_PAR::Vector{FT} = Λ[IΛ_PAR]
    "Wavelength bins for SIF `[nm]`"
    Λ_SIF::Vector{FT} = Λ[IΛ_SIF]
    "Wavelength bins for SIF excitation `[nm]`"
    Λ_SIFE::Vector{FT} = Λ[IΛ_SIFE]
end
