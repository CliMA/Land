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
- [`BetaParameterVcmax`](@ref)
"""
abstract type AbstractBetaParameter end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct to tune G1
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
#     2022-Jun-30: add struct to base on Kleaf
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty struct to indicate to tune G1 or Vcmax based on Kleaf
"""
struct BetaParameterKleaf <: AbstractBetaParameter end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct to base on Ksoil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty struct to indicate to tune G1 or Vcmax based on Ksoil
"""
struct BetaParameterKsoil <: AbstractBetaParameter end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct to base on Pleaf
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty struct to indicate to tune G1 or Vcmax based on Pleaf
"""
struct BetaParameterPleaf <: AbstractBetaParameter end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct to base on Psoil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty struct to indicate to tune G1 or Vcmax based on Psoil
"""
struct BetaParameterPsoil <: AbstractBetaParameter end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct to tune Vcmax
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty struct to indicate to tune Vcmax
"""
struct BetaParameterVcmax <: AbstractBetaParameter end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct to base on Θ (SWC)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty struct to indicate to tune G1 or Vcmax based on Θ (SWC)
"""
struct BetaParameterΘ <: AbstractBetaParameter end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for modular beta function
#     2022-Jun-30: add more types to PARAM_X
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to tune G1 or Vcmax based on leaf hydraulic conductance

# Fields

$(TYPEDFIELDS)

"""
mutable struct BetaFunction
    # parameters that do not change with time
    "Function to turn variables to β tuning factor"
    FUNC::Function
    "Input parameter to base on"
    PARAM_X::Union{BetaParameterKleaf, BetaParameterKsoil, BetaParameterPleaf, BetaParameterPsoil, BetaParameterΘ}
    "Target parameter to tune"
    PARAM_Y::Union{BetaParameterG1, BetaParameterVcmax}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-30: add constructor function
#     2022-Jun-30: fix function definition
#     2022-Jul-07: use BetaParameterKleaf as the default param_x
#
#######################################################################################################################################################################################################
"""

    BetaFunction(f::Function = (func(x) = x); param_x::AbstractBetaParameter = BetaParameterKleaf(), param_y::AbstractBetaParameter = BetaParameterG1())

Construct a `BetaFunction` type beta function, given
- `f` Function
- `param_x` `AbstractBetaParameter` type to indicate which parameter to base on
- `param_y` `AbstractBetaParameter` type to indicate which parameter to tune
"""
BetaFunction(f::Function = (func(x) = x); param_x::AbstractBetaParameter = BetaParameterKleaf(), param_y::AbstractBetaParameter = BetaParameterG1()) = BetaFunction(f, param_x, param_y);


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


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Anderegg model
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for Anderegg stomatal model. The equation used for Anderegg type model is
```math
\\dfrac{∂Θ}{∂E} = \\dfrac{2aP + b}{K}
```
where K is ``\\dfrac{∂E}{∂P}``.

# Fields

$(TYPEDFIELDS)

"""
mutable struct AndereggSM{FT} <: AbstractStomataModel{FT}
    "Quadratic equation parameter `[μmol m⁻² s⁻¹ MPa⁻²]`"
    A::FT
    "Quadratic equation parameter `[μmol m⁻² s⁻¹ MPa⁻¹]`"
    B::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-30: add constructor function
#
#######################################################################################################################################################################################################
"""

    AndereggSM{FT}() where {FT<:AbstractFloat}

Construct a `AndereggSM` type stomatal model
"""
AndereggSM{FT}() where {FT<:AbstractFloat} = AndereggSM{FT}(0.5, 2);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Ball Berry model
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for Ball Berry stomatal model. The equation used for Ball-Berry type model is
```math
gs = g0 + g1 ⋅ RH ⋅ \\dfrac{A}{Cs}
```

# Fields

$(TYPEDFIELDS)

"""
mutable struct BallBerrySM{FT} <: AbstractStomataModel{FT}
    # parameters that do not change with time
    "Minimal stomatal conductance `[mol m⁻² s⁻¹]`"
    G0::FT
    "Slope of conductance-photosynthesis correlation `[-]`"
    G1::FT
    "Beta function to force stomatal response to soil moisture"
    Β::BetaFunction
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-30: add constructor function
#
#######################################################################################################################################################################################################
"""

    BallBerrySM{FT}() where {FT<:AbstractFloat}

Construct a `BallBerrySM` type stomatal model
"""
BallBerrySM{FT}() where {FT<:AbstractFloat} = BallBerrySM{FT}(0.025, 9, BetaFunction());


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Eller model
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty struct for Eller stomatal model. The equation used for Eller type model is
```math
\\dfrac{∂Θ}{∂E} = -\\dfrac{∂K}{∂E} ⋅ \\dfrac{A}{K}
```
where K is ``\\dfrac{∂E}{∂P}``.
"""
struct EllerSM{FT} <: AbstractStomataModel{FT} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Gentine model
#     2022-Jul-07: add field Β to generalize the model
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for Gentine stomatal model. The equation used for Gentine type model is
```math
gs = g0 + g1 ⋅ \\dfrac{k_{leaf}}{k_{max}} ⋅ \\dfrac{A}{Ci}.
```

# Fields

$(TYPEDFIELDS)

"""
mutable struct GentineSM{FT} <: AbstractStomataModel{FT}
    # parameters that do not change with time
    "Minimal stomatal conductance `[mol m⁻² s⁻¹]`"
    G0::FT
    "Slope of conductance-photosynthesis correlation `[-]`"
    G1::FT
    "Beta function to force stomatal response to soil moisture"
    Β::BetaFunction
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-30: add constructor function
#     2022-Jul-07: add field Β to generalize the model (fix the function, param_x, and param_y types)
#
#######################################################################################################################################################################################################
"""

    GentineSM{FT}() where {FT<:AbstractFloat}

Construct a `GentineSM` type stomatal model
"""
GentineSM{FT}() where {FT<:AbstractFloat} = GentineSM{FT}(0.025, 9, BetaFunction(x -> x; param_x = BetaParameterKleaf(), param_y = BetaParameterG1()));


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Leuning model
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for Leuning stomatal model. The equation used for Leuning type model is
```math
gs = g0 + g1 ⋅ \\dfrac{A}{Cs - Γ^{*}} ⋅ \\dfrac{1}{1 + \\dfrac{VPD}{d0}}
```

# Fields

$(TYPEDFIELDS)

"""
mutable struct LeuningSM{FT} <: AbstractStomataModel{FT}
    # parameters that do not change with time
    "Fitting parameter of d/d0 below the fraction, same unit as vpd `[Pa]`"
    D0::FT
    "Minimal stomatal conductance `[mol m⁻² s⁻¹]`"
    G0::FT
    "Slope of conductance-photosynthesis correlation `[-]`"
    G1::FT
    "Beta function to force stomatal response to soil moisture"
    Β::BetaFunction
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-30: add constructor function
#
#######################################################################################################################################################################################################
"""

    LeuningSM{FT}() where {FT<:AbstractFloat}

Construct a `LeuningSM` type stomatal model
"""
LeuningSM{FT}() where {FT<:AbstractFloat} = LeuningSM{FT}(3000, 0.025, 8, BetaFunction());


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Medlyn model
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for Medlyn stomatal model. The equation used for Medlyn type model is
```math
gs = g0 + 1.6 ⋅ \\left( 1 + \\dfrac{g1}{\\sqrt{VPD}} \\right) ⋅ \\dfrac{A}{Ca}
```

# Fields

$(TYPEDFIELDS)

"""
mutable struct MedlynSM{FT} <: AbstractStomataModel{FT}
    # parameters that do not change with time
    "Minimal stomatal conductance `[mol m⁻² s⁻¹]`"
    G0::FT
    "Slope of conductance-photosynthesis correlation `[sqrt(Pa)]`"
    G1::FT
    "Beta function to force stomatal response to soil moisture"
    Β::BetaFunction
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-30: add constructor function
#
#######################################################################################################################################################################################################
"""

    MedlynSM{FT}() where {FT<:AbstractFloat}

Construct a `MedlynSM` type stomatal model
"""
MedlynSM{FT}() where {FT<:AbstractFloat} = MedlynSM{FT}(0.025, 125, BetaFunction());


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Sperry model
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty struct for Sperry stomatal model. The equation used for Sperry type model is
```math
\\dfrac{∂Θ}{∂E} = -\\dfrac{∂K}{∂E} ⋅ \\dfrac{A_{max}}{K_{max}}
```
where K is ``\\dfrac{∂E}{∂P}``.
"""
struct SperrySM{FT} <: AbstractStomataModel{FT} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Wang model
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty struct for Wang stomatal model. The equation used for Wang type model is
```math
\\dfrac{∂Θ}{∂E} = \\dfrac{A}{E_{crit} - E}
```
"""
struct WangSM{FT} <: AbstractStomataModel{FT} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Wang model modified from Anderegg model
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty struct for a new Wang stomatal model modified from Anderegg model. The equation used for new Wang2SM type model is
```math
\\dfrac{∂Θ}{∂E} = \\dfrac{aAP}{K}
```
where K is ``\\dfrac{∂E}{∂P}``.

# Fields

$(TYPEDFIELDS)

"""
mutable struct Wang2SM{FT} <: AbstractStomataModel{FT}
    "Quadratic equation parameter `[MPa⁻²]`"
    A::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-30: add constructor function
#
#######################################################################################################################################################################################################
"""

    Wang2SM{FT}() where {FT<:AbstractFloat}

Construct a `Wang2SM` type stomatal model
"""
Wang2SM{FT}() where {FT<:AbstractFloat} = Wang2SM{FT}(0.1);
