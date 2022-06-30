#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-30: add abstract type for which parameter to tune
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of `AbstractBetaParameter`:
- [`BetaParameterG1`](@ref)
- [`BetaParameterG1`](@ref)
"""
abstract type AbstractBetaParameter end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct to tune stomatal response based on G1
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty struct to indicate to tune G1
"""
struct BetaParameterG1 <: AbstractBetaParameter end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct to tune stomatal response based on Vcmax
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty struct to indicate to tune Vcmax
"""
struct BetaParameterVcmax <: AbstractBetaParameter end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-30: add abstract type for stomatal models
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of `AbstractBetaFunction`:
- [`BetaFunctionLeafK`](@ref)
- [`BetaFunctionSoilK`](@ref)
"""
abstract type AbstractBetaFunction end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct to tune stomatal response based on leaf hydraulic conductance
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to tune G1 or Vcmax based on leaf hydraulic conductance

# Fields

$(TYPEDFIELDS)

"""
mutable struct BetaFunctionLeafK <: AbstractBetaFunction
    # parameters that do not change with time
    "Function to turn variables to β tuning factor"
    FUNC::Function
    "Tuning Vcmax (true) instead of G1"
    PARAM::Union{BetaParameterG1, BetaParameterVcmax}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-30: add constructor function
#
#######################################################################################################################################################################################################
"""

    BetaFunctionLeafK(param::AbstractBetaParameter = BetaParameterG1())

Construct a `BetaFunctionLeafK` type beta function, given
    - `param` `AbstractBetaParameter` type to indicate which parameter to tune
"""
BetaFunctionLeafK(param::AbstractBetaParameter = BetaParameterG1()) = (
    @inline _f(x) = x;

    return BetaFunctionLeafK(_f, param)
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct to tune stomatal response based on soil hydraulic conductance
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to tune G1 or Vcmax based on soil hydraulic conductance

# Fields

$(TYPEDFIELDS)

"""
mutable struct BetaFunctionSoilK <: AbstractBetaFunction
    # parameters that do not change with time
    "Function to turn variables to β tuning factor"
    FUNC::Function
    "Tuning Vcmax (true) instead of G1"
    PARAM::Union{BetaParameterG1, BetaParameterVcmax}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-30: add constructor function
#
#######################################################################################################################################################################################################
"""

    BetaFunctionSoilK(param::AbstractBetaParameter = BetaParameterG1())

Construct a `BetaFunctionSoilK` type beta function, given
- `param` `AbstractBetaParameter` type to indicate which parameter to tune
"""
BetaFunctionSoilK(param::AbstractBetaParameter = BetaParameterG1()) = (
    @inline _f(x) = x;

    return BetaFunctionSoilK(_f, param)
);


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-30: add abstract type for stomatal models
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractStomataModel:
"""
abstract type AbstractStomataModel{FT<:AbstractFloat} end



mutable struct BallBerrySM{FT} <: AbstractStomataModel{FT}
    # parameters that do not change with time
    "Minimal stomatal conductance `[mol m⁻² s⁻¹]`"
    G0::FT
    "Slope of conductance-photosynthesis correlation `[-]`"
    G1::FT
    "Beta function to force stomatal response tp soil moisture"
    Β::AbstractBetaFunction
end
