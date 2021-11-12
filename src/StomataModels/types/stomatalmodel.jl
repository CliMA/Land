###############################################################################
#
# Stomatal model scheme types
#
###############################################################################
#= EmpiricalStomatalModel type tree
AbstractStomatalModel
---> EmpiricalStomatalModel
    ---> ESMBallBerry    # use RH
    ---> ESMGentine      # use k_leaf
    ---> ESMLeuning      # Use VPD
    ---> ESMMedlyn       # use VPD
---> OptimizationStomatalModel
    ---> OSMEller        # Eller 2018 Model
    ---> OSMSperry       # Sperry 2017 model
    ---> OSMWang         # Wang 2020 model
    ---> OSMWAP          # Wolf-Anderegg-Pacala model
    ---> OSMWAP          # Modified Wolf-Anderegg-Pacala model
=#
"""
    type AbstractStomatalModel

Hierarchy of the `AbstractStomatalModel`:
- [`EmpiricalStomatalModel`](@ref)
- [`OptimizationStomatalModel`](@ref)
"""
abstract type AbstractStomatalModel{FT<:AbstractFloat} end




"""
    type EmpiricalStomatalModel

Hierarchy of the `EmpiricalStomatalModel`:
- [`ESMBallBerry`](@ref)
- [`ESMGentine`](@ref)
- [`ESMLeuning`](@ref)
- [`ESMMedlyn`](@ref)
"""
abstract type EmpiricalStomatalModel{FT} <: AbstractStomatalModel{FT} end




"""
    struct ESMBallBerry{FT}

An empirical model parameter set type for Ball-Berry type model.
The equation used for Ball-Berry type model is
```math
gs = g0 + g1 ⋅ RH ⋅ \\dfrac{A}{Cs}
```

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct ESMBallBerry{FT} <: EmpiricalStomatalModel{FT}
    "minimal stomatal conductance g0 `[mol m⁻² s⁻¹]`"
    g0::FT = FT(0.025)
    "slope of conductance-photosynthesis correlation `[unitless]`"
    g1::FT = FT(9.0  )
end




"""
    struct ESMGentine{FT}

An empirical model parameter set type for Gentine type model.
The equation used for Gentine type model is
```math
gs = g0 + g1 ⋅ \\dfrac{k_{leaf}}{k_{max}} ⋅ \\dfrac{A}{Ca}.
```
Note it that the Gentine model does not require for a `β` function to tune the
    soil drought response, but the use of `k_leaf` also does not permit
    post-drought stomatal response unless `k_leaf` can be recovered.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct ESMGentine{FT} <: EmpiricalStomatalModel{FT}
    "minimal stomatal conductance g0 `[mol m⁻² s⁻¹]`"
    g0::FT = FT(0.025)
    "slope of conductance-photosynthesis correlation `[unitless]`"
    g1::FT = FT(9.0  )
end




"""
    struct ESMLeuning{FT}

An empirical model parameter set type for Leuning type model.
The equation used for Leuning type model is
```math
gs = g0 + g1 ⋅ \\dfrac{A}{Cs - Γ^{*}} ⋅ \\dfrac{1}{1 + \\dfrac{VPD}{d0}}
```

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct ESMLeuning{FT} <: EmpiricalStomatalModel{FT}
    "minimal stomatal conductance g0 `[mol m⁻² s⁻¹]`"
    g0::FT = FT(0.025 )
    "slope of conductance-photosynthesis correlation `[unitless]`"
    g1::FT = FT(8.0   )
    "fitting parameter of d/d0 below the fraction, same unit as vpd `[Pa]`"
    d0::FT = FT(3000.0)
end




"""
    struct ESMMedlyn{FT}

An empirical model parameter set type for Medlyn type model.
The equation used in Medlyn type model is
```math
gs = g0 + 1.6 ⋅ \\left( 1 + \\dfrac{g1}{\\sqrt{VPD}} \\right) ⋅ \\dfrac{A}{Ca}
```

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct ESMMedlyn{FT} <: EmpiricalStomatalModel{FT}
    "minimal stomatal conductance g0 `[mol m⁻² s⁻¹]`"
    g0::FT = FT(0.025)
    "slope of conductance-photosynthesis correlation `[Pa⁽⁵⁾]`"
    g1::FT = FT(125.0)
end




"""
    type OptimizationStomatalModel

Hierarchy of the `OptimizationStomatalModel`:
- [`OSMEller`](@ref)
- [`OSMSperry`](@ref)
- [`OSMWang`](@ref)
- [`OSMWAP`](@ref)
- [`OSMWAPMod`](@ref)
"""
abstract type OptimizationStomatalModel{FT} <: AbstractStomatalModel{FT} end




"""
    struct OSMEller

An optimization model parameter set type for Eller model.
The equation used for Eller model is
```math
\\dfrac{∂Θ}{∂E} = -\\dfrac{∂K}{∂E} ⋅ \\dfrac{A}{K}
```
where K is ``\\dfrac{∂E}{∂P}``.

# Fields
$(TYPEDFIELDS)
"""
struct OSMEller{FT} <: OptimizationStomatalModel{FT} end




"""
    struct OSMSperry

An optimization model parameter set type for Sperry model.
The equation used for Sperry model is
```math
\\dfrac{∂Θ}{∂E} = -\\dfrac{∂K}{∂E} ⋅ \\dfrac{A_{max}}{K_{max}}
```
where K is ``\\dfrac{∂E}{∂P}``.

# Fields
$(TYPEDFIELDS)
"""
struct OSMSperry{FT} <: OptimizationStomatalModel{FT} end




"""
    struct OSMWang

An optimization model parameter set type for Eller type model.
The equation used for Wang model is
```math
\\dfrac{∂Θ}{∂E} = \\dfrac{A}{E_{crit} - E}
```

# Fields
$(TYPEDFIELDS)
"""
struct OSMWang{FT} <: OptimizationStomatalModel{FT} end




"""
    struct OSMWAP{FT}

An optimization model parameter set type for Wolf-Anderegg-Pacala type model.
The equation used for Wolf-Anderegg-Pacala model is
```math
\\dfrac{∂Θ}{∂E} = \\dfrac{2aP + b}{K}
```
where K is ∂P/∂E.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct OSMWAP{FT} <: OptimizationStomatalModel{FT}
    "Quadratic equation parameter `[μmol m⁻² s⁻¹ MPa⁻²]`"
    a::FT = FT(0.5)
    "Quadratic equation parameter `[μmol m⁻² s⁻¹ MPa⁻¹]`"
    b::FT = FT(2.0)
end




"""
    struct OSMWAPMod{FT}

An optimization model parameter set type for Wolf-Anderegg-Pacala type model,
    modified by adding a photosynthesis component while set b and c = 0.
The equation used for modified Wolf-Anderegg-Pacala model is
```math
\\dfrac{∂Θ}{∂E} = \\dfrac{aAP}{K}
```
where P is absolute value of leaf xylem pressure.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct OSMWAPMod{FT} <: OptimizationStomatalModel{FT}
    "Quadratic equation parameter `[mol mol⁻¹ MPa⁻¹]`"
    a::FT = FT(0.1)
end
