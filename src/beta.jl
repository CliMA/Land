#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-30: migrate function from older version
#
#######################################################################################################################################################################################################
"""
This function returns the beta function to force stomatal response to soil moisture. Supported methods are

$(METHODLIST)

"""
function β_factor end



β_factor(f::Function, param::BetaParameterKleaf, hs::LeafHydraulics, p_leaf::FT) where {FT<:AbstractFloat} = f(relative_hydraulic_conductance(hs.VC, p_leaf));



β_factor(f::Function, param::BetaParameterKsoil, vc::AbstractSoilVC, ψ::Bool, ψ_soil::FT) where {FT<:AbstractFloat} = f(relative_hydraulic_conductance(vc, ψ, ψ_soil));



β_factor(f::Function, param::Union{BetaParameterPleaf, BetaParameterPsoil, BetaParameterΘ}, x::FT) where {FT<:AbstractFloat} = f(x);
