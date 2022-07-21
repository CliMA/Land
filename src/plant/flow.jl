#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-27: add abstract type for steady and non-steady flow
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractFlowProfile:
- [`NonSteadyStateFlow`](@ref)
- [`SteadyStateFlow`](@ref)

"""
abstract type AbstractFlowProfile{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-27: add non-steady state flow struct
#     2022-May-31: remove unneeded variables from struct
#     2022-May-31: add field: f_sum
#     2022-Jul-19: add dimension control to struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains stem hydraulic system flow rates at non-steady state

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct NonSteadyStateFlow{FT<:AbstractFloat} <: AbstractFlowProfile{FT}
    # Dimensions
    "Dimension of capaciatance elements"
    DIM_CAPACITY::Int = 1

    # Diagnostic variables
    "Flow rate into the organ `[mol s⁻¹]` (for root and stem) or `[mol m⁻² s⁻¹]` (for leaf)"
    f_in::FT = 0
    "Flow rate out of the organ `[mol s⁻¹]` (for root and stem) or `[mol m⁻² s⁻¹]` (for leaf)"
    f_out::FT = 0

    # Cache variables
    "Vector of xylem water flow `[mol s⁻¹]` (for root and stem) or `[mol m⁻² s⁻¹]` (for leaf)"
    _f_element::Vector{FT} = zeros(FT, DIM_CAPACITY)
    "Vector of buffer water flow `[mol s⁻¹]` (for root and stem) or `[mol m⁻² s⁻¹]` (for leaf)"
    _f_buffer::Vector{FT} = zeros(FT, DIM_CAPACITY)
    "Vector of sum buffer water flow `[mol s⁻¹]` (for root and stem) or `[mol m⁻² s⁻¹]` (for leaf)"
    _f_sum::Vector{FT} = zeros(FT, DIM_CAPACITY)
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-27: add non-steady state flow struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains stem hydraulic system flow rates at steady state

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SteadyStateFlow{FT<:AbstractFloat} <: AbstractFlowProfile{FT}
    # Diagnostic variables
    "Flow rate through the organ `[mol s⁻¹]` (for root and stem) or `[mol m⁻² s⁻¹]` (for leaf)"
    flow::FT = 0
end
