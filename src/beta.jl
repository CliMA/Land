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


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Junl-01: add method to tune stomatal opening based on relative hydraulic conductance at leaf xylem end
#
#######################################################################################################################################################################################################
"""

    β_factor(f::Function, hs::LeafHydraulics{FT}, p_25::FT) where {FT<:AbstractFloat}

Return the β factor based on relative leaf hydraulic conductance at xylem end, given
- `f` Function to translate relative k to β, for example f(x) = x, f(x) = x², and f(x) = sqrt(x) for x in [0,1]
- `hs` Leaf hydraulics system
- `p_25` Leaf xylem pressure at the end of xylem (corrected to 25 °C)
"""
β_factor(f::Function, hs::LeafHydraulics{FT}, p_25::FT) where {FT<:AbstractFloat} = FT(f(relative_hydraulic_conductance(hs.VC, p_25)));


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Junl-01: add method to tune stomatal opening based on relative hydraulic conductance of the soil
#
#######################################################################################################################################################################################################
"""

    β_factor(f::Function, vc::AbstractSoilVC{FT}, ψ_25::FT) where {FT<:AbstractFloat}

Return the β factor based on relative hydraulic conductance of the soil, given
- `f` Function to translate relative k to β, for example f(x) = x, f(x) = x², and f(x) = sqrt(x) for x in [0,1]
- `vc` Soil vulnerability curve (moisture retention curve)
- `ψ_25` Soil metric potential corrected to 25 °C (note that this function may not be useful for plants with salt stress)
"""
β_factor(f::Function, vc::AbstractSoilVC{FT}, ψ_25::FT) where {FT<:AbstractFloat} = FT(f(relative_hydraulic_conductance(vc, true, ψ_25)));


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Junl-01: add method to tune stomatal opening based on soil potential or leaf pressure
#
#######################################################################################################################################################################################################
"""

    β_factor(f::Function, x_25::FT) where {FT<:AbstractFloat}

Return the β factor based on soil water potential, leaf xylem pressure, or soil water content, given
- `f` Function to translate x to β, for example f(x) = x, f(x) = x², and f(x) = sqrt(x) for x in [0,1]
- `x_25` Leaf xylem pressure corrected to 25 °C, soil water potential corrected to 25 °C (forcing on roots), or soil water content
"""
β_factor(f::Function, x_25::FT) where {FT<:AbstractFloat} = β_factor(f, x_25);



#=
# x for θ
β_factor(f::Function, θ::FT) where {FT<:AbstractFloat} = FT(f(θ));



β_factor(f::Function, param::BetaParameterKleaf, hs::LeafHydraulics{FT}, p_leaf::FT, ψ_soil::FT, θ::FT) where {FT<:AbstractFloat} = β_factor(f, hs, p_leaf)



β_factor(f::Function, param::BetaParameterPleaf, hs::LeafHydraulics{FT}, p_leaf::FT, ψ_soil::FT, θ::FT) where {FT<:AbstractFloat} = β_factor(f, p_leaf);



β_factor(f::Function, param::BetaParameterKsoil, hs::LeafHydraulics{FT}, vc::AbstractSoilVC{FT}, p_leaf::FT, ψ_soil::FT, θ::FT) where {FT<:AbstractFloat} = β_factor(f, vc, ψ_soil);



β_factor(f::Function, param::BetaParameterPsoil, hs::LeafHydraulics{FT}, vc::AbstractSoilVC{FT}, p_leaf::FT, ψ_soil::FT, θ::FT) where {FT<:AbstractFloat} = β_factor(f, ψ_soil);



β_factor(f::Function, param::BetaParameterΘ, hs::LeafHydraulics{FT}, vc::AbstractSoilVC{FT}, p_leaf::FT, ψ_soil::FT, θ::FT) where {FT<:AbstractFloat} = β_factor(f, θ);
=#
