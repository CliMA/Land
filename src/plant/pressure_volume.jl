#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Apr-20: add abstract type for pressure volume curve
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierachy of AbstractPVCurve:
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
mutable struct LinearPVCurve{FT<:AbstractFloat} <: AbstractPVCurve{FT}
    "Conductance for refilling (relative to maximum) `[MPa⁻¹ s⁻¹]`"
    K_REFILL::FT
    "Slope of the linear PV curve (relative to maximum) `[MPa⁻¹]`"
    SLOPE::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-may-24: add default constructor
#
#######################################################################################################################################################################################################
"""

    LinearPVCurve{FT}() where {FT<:AbstractFloat}

Constructor for LinearPVCurve

---
# Examples
```julia
pvc = LinearPVCurve{Float64}();
```
"""
LinearPVCurve{FT}() where {FT<:AbstractFloat} = LinearPVCurve{FT}(0.0001, 0.2);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-24: add segmented PV curve
#     2022-May-25: fix floating number type control
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains information for segmented PV curve

# Fields

$(TYPEDFIELDS)

"""
mutable struct SegmentedPVCurve{FT} <: AbstractPVCurve{FT}
    "n_o / maximum V `[mol m⁻³]`"
    C_ALL::FT
    "Conductance for refilling (relative to maximum) `[MPa⁻¹ s⁻¹]`"
    K_REFILL::FT
    "Apoplastic water content relative to maximum water volume"
    RWC_APO::FT
    "Relative water content at turgor loss point"
    RWC_TLP::FT
    "Bulk modulus of elasticity `[MPa]`"
    Ε_BULK::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-may-24: add default constructor
#
#######################################################################################################################################################################################################
"""

    SegmentedPVCurve{FT}() where {FT<:AbstractFloat}

Constructor for SegmentedPVCurve

---
# Examples
```julia
pvc = SegmentedPVCurve{Float64}();
```
"""
SegmentedPVCurve{FT}() where {FT<:AbstractFloat} = SegmentedPVCurve{FT}(300.0, 0.0001, 0.2, 0.8, 20.0);
