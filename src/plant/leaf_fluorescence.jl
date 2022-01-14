"""
$(TYPEDEF)

Hierachy of the `AbstractFluorescenceModel`
- [``](@ref)
"""
abstract type AbstractFluorescenceModel{FT<:AbstractFloat} end


"""
$(TYPEDEF)

Structure that stores van der Tol et al. (2013) fluorescence model parameters.

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
FluorescenceVDT(FT) = VanDerTolFluorescenceModel{FT}(2.48, 2.83, 0.114)
FluorescenceVDTDrought(FT) = VanDerTolFluorescenceModel{FT}(5.01, 1.93, 10)
