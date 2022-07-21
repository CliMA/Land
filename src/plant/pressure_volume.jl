#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Apr-20: add abstract type for pressure volume curve
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractPVCurve:
- [`LinearPVCurve`](@ref)
- [`SegmentedPVCurve`](@ref)

"""
abstract type AbstractPVCurve{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Apr-20: add linear PV curve
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains information for linear PV curve

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LinearPVCurve{FT<:AbstractFloat} <: AbstractPVCurve{FT}
    # General model information
    "Conductance for refilling (relative to maximum) `[MPa⁻¹ s⁻¹]`"
    K_REFILL::FT = 1e4
    "Slope of the linear PV curve (relative to maximum) `[MPa⁻¹]`"
    SLOPE::FT = 0.2
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-24: add segmented PV curve
#     2022-Jul-20: rename Ε_BULK (greek) to ϵ_BULK
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains information for segmented PV curve

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SegmentedPVCurve{FT<:AbstractFloat} <: AbstractPVCurve{FT}
    # General model information
    "n_o / maximum V `[mol m⁻³]`"
    C_ALL::FT = 300
    "Conductance for refilling (relative to maximum) `[MPa⁻¹ s⁻¹]`"
    K_REFILL::FT = 1e-4
    "Apoplastic water content relative to maximum water volume"
    RWC_APO::FT = 0.2
    "Relative water content at turgor loss point"
    RWC_TLP::FT = 0.8
    "Bulk modulus of elasticity `[MPa]`"
    ϵ_BULK::FT = 20
end
