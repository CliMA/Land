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
#     2022-Jul-08: use @kwdef for the constructor
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
Base.@kwdef mutable struct AndereggSM{FT} <: AbstractStomataModel{FT}
    "Quadratic equation parameter `[μmol m⁻² s⁻¹ MPa⁻²]`"
    A::FT = 0.5
    "Quadratic equation parameter `[μmol m⁻² s⁻¹ MPa⁻¹]`"
    B::FT = 2
    "Slope constant `[mol² m⁻² s⁻¹ μmol⁻¹]`"
    K::FT = 1e-7
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Ball Berry model
#     2022-Jul-07: add time constant
#     2022-Jul-08: use @kwdef for the constructor
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
Base.@kwdef mutable struct BallBerrySM{FT} <: AbstractStomataModel{FT}
    # parameters that do not change with time
    "Minimal stomatal conductance `[mol m⁻² s⁻¹]`"
    G0::FT = 0.025
    "Slope of conductance-photosynthesis correlation `[-]`"
    G1::FT = 9
    "Beta function to force stomatal response to soil moisture"
    Β::BetaFunction = BetaFunction()
    "Time constant for the prognostic stomatal conductance `[s]`"
    Τ::FT = 600
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Eller model
#     2022-Jul-08: use @kwdef for the constructor
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
Base.@kwdef mutable struct EllerSM{FT} <: AbstractStomataModel{FT}
    "Slope constant `[mol² m⁻² s⁻¹ μmol⁻¹]`"
    K::FT = 1e-7
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Gentine model
#     2022-Jul-07: add field Β to generalize the model
#     2022-Jul-07: add time constant
#     2022-Jul-08: use @kwdef for the constructor
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
Base.@kwdef mutable struct GentineSM{FT} <: AbstractStomataModel{FT}
    # parameters that do not change with time
    "Minimal stomatal conductance `[mol m⁻² s⁻¹]`"
    G0::FT = 0.025
    "Slope of conductance-photosynthesis correlation `[-]`"
    G1::FT = 9
    "Beta function to force stomatal response to soil moisture"
    Β::BetaFunction = BetaFunction(x -> x; param_x = BetaParameterKleaf(), param_y = BetaParameterG1())
    "Time constant for the prognostic stomatal conductance `[s]`"
    Τ::FT = 600
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Leuning model
#     2022-Jul-07: add time constant
#     2022-Jul-08: use @kwdef for the constructor
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
Base.@kwdef mutable struct LeuningSM{FT} <: AbstractStomataModel{FT}
    # parameters that do not change with time
    "Fitting parameter of d/d0 below the fraction, same unit as vpd `[Pa]`"
    D0::FT = 3000
    "Minimal stomatal conductance `[mol m⁻² s⁻¹]`"
    G0::FT = 0.025
    "Slope of conductance-photosynthesis correlation `[-]`"
    G1::FT = 8
    "Beta function to force stomatal response to soil moisture"
    Β::BetaFunction = BetaFunction()
    "Time constant for the prognostic stomatal conductance `[s]`"
    Τ::FT = 600
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Medlyn model
#     2022-Jul-07: add time constant
#     2022-Jul-08: use @kwdef for the constructor
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
Base.@kwdef mutable struct MedlynSM{FT} <: AbstractStomataModel{FT}
    # parameters that do not change with time
    "Minimal stomatal conductance `[mol m⁻² s⁻¹]`"
    G0::FT = 0.025
    "Slope of conductance-photosynthesis correlation `[sqrt(Pa)]`"
    G1::FT = 125
    "Beta function to force stomatal response to soil moisture"
    Β::BetaFunction = BetaFunction()
    "Time constant for the prognostic stomatal conductance `[s]`"
    Τ::FT = 600
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Sperry model
#     2022-Jul-08: use @kwdef for the constructor
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
Base.@kwdef mutable struct SperrySM{FT} <: AbstractStomataModel{FT}
    "Slope constant `[mol² m⁻² s⁻¹ μmol⁻¹]`"
    K::FT = 1e-7
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Wang model
#     2022-Jul-08: use @kwdef for the constructor
#     2022-Jul-11: add field fitness factor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Empty struct for Wang stomatal model. The equation used for Wang type model is
```math
\\dfrac{∂Θ}{∂E} = \\dfrac{A}{E_{crit} - E}
```
"""
Base.@kwdef mutable struct WangSM{FT} <: AbstractStomataModel{FT}
    # parameters that do not change with time
    "Fitness factor"
    F_FITNESS::FT = 0.1
    "Slope constant `[mol² m⁻² s⁻¹ μmol⁻¹]`"
    K::FT = 1e-7

    # diagnostic variables that change with time
    "Ratio that leaf area is exposed to external sources/sinks (not other leaves, e.g., 2/LAI for canopy on average, used for nocturnal transpiration)"
    f_view::FT = 2
    "Memory PPAR `[μmol m⁻² s⁻¹]`"
    ppar_mem::FT = 100
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-30: add struct for Wang model modified from Anderegg model
#     2022-Jul-08: use @kwdef for the constructor
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
Base.@kwdef mutable struct Wang2SM{FT} <: AbstractStomataModel{FT}
    "Quadratic equation parameter `[MPa⁻²]`"
    A::FT = 0.1
    "Slope constant `[mol² m⁻² s⁻¹ μmol⁻¹]`"
    K::FT = 1e-7
end
