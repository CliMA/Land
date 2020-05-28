#= AbstractEmpiricalStomatalModel type tree
AbstractStomatalModel
---> AbstractEmpiricalStomatalModel
    ---> ESMBallBerry    # use RH
    ---> ESMGentine
    ---> ESMLeuning
    ---> ESMMedlyn       # use VPD
=#

abstract type AbstractStomatalModel end

abstract type AbstractEmpiricalStomatalModel <: AbstractStomatalModel end




"""
    ESMBallBerry{FT}

An empirical model parameter set type for Ball-Berry type model.
The equation used for Ball-Berry type model is `gs = g0 + g1 * RH * A/Ci`.
Note it that these empirical model often takes Ci or Cc as input, so a conversion from Pa to ppm is required.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ESMBallBerry{FT} <: AbstractEmpiricalStomatalModel
    "minimal stomatal conductance g0 `[mol m⁻² s⁻¹]`, should this one be temperature-dependent?"
    g0::FT = FT(0.025)
    "slope of conductance-photosynthesis correlation `[unitless]`"
    g1::FT = FT(9.0  )
end




"""
    ESMGentine{FT}

An empirical model parameter set type for Gentine type model.
The equation used for Gentine type model is `gs = g0 + g1 * k_leaf/k_max * A/Ci`.
Note it that these empirical model often takes Ci or Cc as input, so a conversion from Pa to ppm is required.
Note it that the Gentine model does not require for a `β` function to tune the soil drought response, but the use of `k_leaf` also does not permit post-drought stomatal response unless `k_leaf` can be recovered. 

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ESMGentine{FT} <: AbstractEmpiricalStomatalModel
    "minimal stomatal conductance g0 `[mol m⁻² s⁻¹]`, should this one be temperature-dependent?"
    g0::FT = FT(0.025)
    "slope of conductance-photosynthesis correlation `[unitless]`"
    g1::FT = FT(9.0  )
end




"""
    ESMLeuning{FT}

An empirical model parameter set type for Leuning type model.
The equation used for Leuning type model is `gs = g0 + g1 * A/(Ci-Γ_star) * 1/(1+VPD/d0)`.
Note it that these empirical model often takes Ci or Cc as input, so a conversion from Pa to ppm is required.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ESMLeuning{FT} <: AbstractEmpiricalStomatalModel
    "minimal stomatal conductance g0 `[mol m⁻² s⁻¹]`, should this one be temperature-dependent?"
    g0::FT = FT(0.025 )
    "slope of conductance-photosynthesis correlation `[unitless]`"
    g1::FT = FT(8.0   )
    "fitting parameter of d/d0 below the fraction, same unit as vpd `[Pa]`"
    d0::FT = FT(3000.0)
end




"""
    ESMMedlyn{FT}

An empirical model parameter set type for Medlyn type model.
The equation used in Medlyn type model is `gs = g0 + (1+g1/sqrt(VPD)) * A/Ci`.
Note it that these empirical model often takes Ci or Cc as input, so a conversion from Pa to ppm is required.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ESMMedlyn{FT} <: AbstractEmpiricalStomatalModel
    "minimal stomatal conductance g0 `[mol m⁻² s⁻¹]`, should this one be temperature-dependent?"
    g0::FT = FT(0.025)
    "slope of conductance-photosynthesis correlation `[Pa⁽⁵⁾]`"
    g1::FT = FT(125.0)
end
