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
# Changes to this type
# General
#     2022-May-27: add non-steady state flow struct
#     2022-May-31: remove unneeded variables from struct
#     2022-May-31: add f_sum field
#     2022-Jul-18: use kwdef for the constructor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains stem hydraulic system flow rates at non-steady state

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct NonSteadyStateFlow{FT<:AbstractFloat} <: AbstractFlowProfile{FT}
    # parameters that do not change with time
    "Number of capaciatance elements"
    N::Int = 1

    # dignostic variables that change with time
    "Vector of buffer water flow `[mol m⁻²]`"
    f_buffer::Vector{FT} = zeros(FT, N)
    "Vector of xylem water flow `[mol m⁻²]`"
    f_element::Vector{FT} = zeros(FT, N)
    "Flow rate in `[mol s⁻¹]` or `[mol m⁻² s⁻¹]` (for leaf)"
    f_in::FT = 0
    "Flow rate out `[mol s⁻¹]` or `[mol m⁻² s⁻¹]` (for leaf)"
    f_out::FT = 0
    "Vector of sum buffer water flow `[mol m⁻²]`"
    f_sum::Vector{FT} = zeros(FT, N)
end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-27: add non-steady state flow struct
#     2022-Jul-18: use kwdef for the constructor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains stem hydraulic system flow rates at steady state

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SteadyStateFlow{FT<:AbstractFloat} <: AbstractFlowProfile{FT}
    "Flow rate `[mol s⁻¹]` or `[mol m⁻² s⁻¹]` (for leaf)"
    flow::FT = 0
end
