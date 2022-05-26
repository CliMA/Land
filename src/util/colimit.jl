#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jan-18: add abstract colimitation type
#     2022-Jan-25: add documentation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of `AbstractColimit`
- [`MinimumColimit`](@ref)
- [`QuadraticColimit`](@ref)
- [`SerialColimit`](@ref)
"""
abstract type AbstractColimit{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-14: add minimum colimitation struct
#     2022-Jan-25: add documentation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty structure to indicate minimum colimitation.

---
# Examples
```julia
col = MinimumColimit{Float64}();
```
"""
struct MinimumColimit{FT<:AbstractFloat} <:AbstractColimit{FT} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-14: add quadratic colimitation struct
#     2022-Jan-25: add documentation
#     2022-Mar-01: add documentation to field
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to indicate quadratic colimitation.

# Fields

$(TYPEDFIELDS)

---
# Examples
```julia
col = QuadraticColimit{Float64}(0.98);
```
"""
mutable struct QuadraticColimit{FT<:AbstractFloat} <: AbstractColimit{FT}
    "Curvature factor"
    CURVATURE::FT
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Feb-28: add serial colimitation struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty structure to indicate serial colimitation.

---
# Examples
```julia
col = SerialColimit{Float64}();
```
"""
struct SerialColimit{FT<:AbstractFloat} <:AbstractColimit{FT} end
