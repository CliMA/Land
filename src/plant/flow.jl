#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-27: add abstract type for steady and non-steady flow
#     2022-May-31: add documentation
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
#     2022-May-31: add documentation
#     2022-May-31: add f_sum field
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains stem hydraulic system flow rates at non-steady state

# Fields

$(TYPEDFIELDS)

"""
mutable struct NonSteadyStateFlow{FT} <: AbstractFlowProfile{FT}
    "Vector of buffer water flow `[mol m⁻²]`"
    f_buffer::Vector{FT}
    "Vector of xylem water flow `[mol m⁻²]`"
    f_element::Vector{FT}
    "Flow rate in `[mol s⁻¹]` or `[mol m⁻² s⁻¹]` (for leaf)"
    f_in::FT
    "Flow rate out `[mol s⁻¹]` or `[mol m⁻² s⁻¹]` (for leaf)"
    f_out::FT
    "Vector of sum buffer water flow `[mol m⁻²]`"
    f_sum::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-27: add non-steady state flow struct
#     2022-May-31: remove unneeded variables from struct
#     2022-May-31: add documentation
#     2022-May-31: add f_sum field
#
#######################################################################################################################################################################################################
"""

    NonSteadyStateFlow{FT}(N::Int, isleaf::Bool = true) where {FT<:AbstractFloat}

Construct a non-steady state flow struct, given
- `N` Number of buffer rates from capaciatance (always 1 for Leaf)
- `isleaf` Bool to indicate if the organ is a leaf
"""
NonSteadyStateFlow{FT}(N::Int, isleaf::Bool = true) where {FT<:AbstractFloat} = (
    if isleaf
        _zeros = zeros(FT,1);
    else
        _zeros = zeros(FT,N);
    end;

    return NonSteadyStateFlow{FT}(
                _zeros,     # f_buffer
                _zeros,     # f_element
                0,          # f_in
                0,          # f_out
                _zeros      # f_sum
    )
);


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-27: add non-steady state flow struct
#     2022-May-31: add documentation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains stem hydraulic system flow rates at steady state

# Fields

$(TYPEDFIELDS)

"""
mutable struct SteadyStateFlow{FT} <: AbstractFlowProfile{FT}
    "Flow rate `[mol s⁻¹]` or `[mol m⁻² s⁻¹]` (for leaf)"
    flow::FT
end
