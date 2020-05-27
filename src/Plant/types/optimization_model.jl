#= AbstractOptimizationStomatalModel type tree
AbstractStomatalModel
---> AbstractOptimizationStomatalModel
    ---> OSMEller     # Eller 2018 Model
    ---> OSMSperry    # Sperry 2017 model
    ---> OSMWang      # Wang 2020 model
    ---> OSMWAP       # Wolf-Anderegg-Pacala model
=#

abstract type AbstractOptimizationStomatalModel <: AbstractStomatalModel end

struct OSMEller  <: AbstractOptimizationStomatalModel end
struct OSMSperry <: AbstractOptimizationStomatalModel end
struct OSMWang   <: AbstractOptimizationStomatalModel end




"""
    OSMWAP{FT}

An optimization model parameter set type for Wolf-Anderegg-Pacala type model.
The equation used for Wolf-Anderegg-Pacala model is `Θ = a*P^2 + b*P +c`, where P is absolute value of leaf xylem pressure.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct OSMWAP{FT}
    "Quadratic equation parameter `[μmol m⁻² s⁻¹ MPa⁻²]`"
    a::FT = FT(0.5)
    "Quadratic equation parameter `[μmol m⁻² s⁻¹ MPa⁻¹]`"
    b::FT = FT(2.0)
end




"""
    OSMWAPMod{FT}

An optimization model parameter set type for Wolf-Anderegg-Pacala type model, modified by adding a photosynthesis component while set b and c = 0.
The equation used for modified Wolf-Anderegg-Pacala model is `∂Θ/∂E = a * P`, where P is absolute value of leaf xylem pressure.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct OSMWAPMod{FT}
    "Quadratic equation parameter `[mol mol⁻¹ MPa⁻¹]`"
    a::FT = FT(0.1)
end
