# Provides a framework on how to translate between LAI coordinates and physical Z-dimensions
#=

using Distributions:Beta,cdf,pdf
using Parameters

abstract type AbstractVerticalCanopyStructure end


"The beta distribution probability for LAI as function of z/h"
Base.@kwdef struct BetaCanopyStructure{FT} <: AbstractVerticalCanopyStructure
    p::FT = 2.5
    q::FT = 2.5
    B = Beta(p,q)
end

# Bonan ML canopy distribution (Beta function):
function convert_z_to_dLAI!(mod::BetaCanopyStructure, canopy, z)
    @unpack LAI, nlayers, h = canopy
    # LAI = Leaf Area Index
    # h = Canopy height
    LAI * cdf.(mod.B,z/h)
end # functi
=#
