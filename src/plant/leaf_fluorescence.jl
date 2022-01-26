#######################################################################################################################################################################################################
#
# Changes to the type
# General
#     2022-Jan-14: add abstract fluorescence model type
#     2022-Jan-24: fix documentation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierachy of the `AbstractFluorescenceModel`
- [`VanDerTolFluorescenceModel`](@ref)
- [`CytochromeFluorescenceModel`](@ref)
"""
abstract type AbstractFluorescenceModel{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to the type
# General
#     2022-Jan-14: add van der Tol model struct
#     2022-Jan-24: fix documentation
# Sources
#     van der Tol et al. (2014) Models of fluorescence and photosynthesis for interpreting measurements of solar-induced chlorophyll fluorescence
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores van der Tol et al. (2014) fluorescence model parameters.

# Fields

$(TYPEDFIELDS)

"""
struct VanDerTolFluorescenceModel{FT<:AbstractFloat} <: AbstractFluorescenceModel{FT}
    "Fitting parameter K_0"
    K_0::FT
    "Fitting parameter α"
    K_A::FT
    "Fitting parameter β"
    K_B::FT
end


#######################################################################################################################################################################################################
#
# Changes to the constructors
# General
#     2022-Jan-14: migrate from Photosynthesis.jl
# Sources
#     van der Tol et al. (2013) Models of fluorescence and photosynthesis for interpreting measurements of solar-induced chlorophyll fluorescence
#
#######################################################################################################################################################################################################
"""

    VanDerTolFluorescenceModel{FT}(drought::Bool = false) where {FT<:AbstractFloat}

Constructor for `VanDerTolFluorescenceModel` fluorescence model, given
- `drought` If true, use parameters trained from drought stressed plant. Default is `false`.

---
# Examples
```julia
vdt = VanDerTolFluorescenceModel{Float64}();
vdt = VanDerTolFluorescenceModel{Float64}(true);
```
"""
VanDerTolFluorescenceModel{FT}(drought::Bool = false) where {FT<:AbstractFloat} = (
    if drought
        return VanDerTolFluorescenceModel{FT}(5.01, 1.93, 10)
    else
        return VanDerTolFluorescenceModel{FT}(2.48, 2.83, 0.114)
    end
);


#######################################################################################################################################################################################################
#
# Changes to the structure
# General
#     2022-Jan-24: add an empty structure for Cytochrome based fluorescence model
#     2022-Jan-24: fix documentation
# Sources
#     Johnson and Berry (2021) The role of Cytochrome b₆f in the control of steady-state photosynthesis: a conceptual and quantitative model
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores Johnson and Berry (2021) fluorescence model (empty structure).

# Fields

$(TYPEDFIELDS)

"""
struct CytochromeFluorescenceModel{FT<:AbstractFloat} <: AbstractFluorescenceModel{FT} end
