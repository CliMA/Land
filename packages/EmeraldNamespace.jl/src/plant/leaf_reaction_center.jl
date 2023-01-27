#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2022-Jan-14: add van der Tol model struct
#     2022-Feb-07: remove the Hierarchy from abstract fluorescence model
#     2022-Feb-07: move struct definition as a field of VJPReactionCenter
#     2022-Oct-27: make VanDerTolFluorescenceModel struct mutable
# To do
#     TODO: examine why van der Tol et al has the nonstressed parameter set that are so far away from the stressed one
# Sources
#     van der Tol et al. (2013) Models of fluorescence and photosynthesis for interpreting measurements of solar-induced chlorophyll fluorescence
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores van der Tol et al. (2014) fluorescence model parameters.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct VanDerTolFluorescenceModel{FT<:AbstractFloat}
    # General model information
    "Fitting parameter K_0"
    K_0::FT = 5.01
    "Fitting parameter α"
    K_A::FT = 1.93
    "Fitting parameter β"
    K_B::FT = 10
end


""" VanDerTolFluorescenceModel that uses data from all observations """
VDTModelAll(FT) = VanDerTolFluorescenceModel{FT}(K_0 = 2.48, K_A = 2.83, K_B = 0.114)


""" VanDerTolFluorescenceModel that uses data from drought stressed observations """
VDTModelDrought(FT) = VanDerTolFluorescenceModel{FT}(K_0 = 5.01, K_A = 1.93, K_B = 10);


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jan-24: abstractize the reaction center
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Abstract type for reaction center

Hierarchy of the `AbstractReactionCenter`
- [`VJPReactionCenter`](@ref)
- [`CytochromeReactionCenter`](@ref)

"""
abstract type AbstractReactionCenter{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-15: isolate the reaction center from Leaf in Photosynthesis.jl
#     2022-Feb-07: add fluorescence model as a field
#     2022-Jul-20: rename quite some cache variables
#     2023-Jan-27: set default K_D to 0.95 instead of 0.85
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores reaction center information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct VJPReactionCenter{FT<:AbstractFloat} <:AbstractReactionCenter{FT}
    # Constant coefficients
    "Fraction of absorbed light used by PSII ETR"
    F_PSII::FT = 0.5
    "Rate constant for thermal dissipation"
    K_D::FT = 0.95
    "Rate constant for fluorescence"
    K_F::FT = 0.05
    "Maximal rate constant for photochemistry"
    K_P_MAX::FT = 4
    "max PSII yield (k_npq_r=0, all RC open)"
    Φ_PSII_MAX::FT = K_P_MAX / (K_D + K_F + K_P_MAX)

    # Embedded structures
    "Fluorescence model"
    FLM::VanDerTolFluorescenceModel{FT} = VanDerTolFluorescenceModel{FT}()

    # Prognostic variables
    "Sustained NPQ rate constant (for seasonal changes, default is zero)"
    k_npq_sus::FT = 0

    # Diagnostic variables
    "Fluorescence yield"
    ϕ_f::FT = 0
    "Photochemical yield"
    ϕ_p::FT = 0

    # Cache variables
    "Dark adapted yield (`Kp=0`)"
    _f_m::FT = 0
    "Light adapted yield (`Kp=0`)"
    _f_m′::FT = 0
    "Dark-adapted fluorescence yield (`Kp=max`)"
    _f_o::FT = 0
    "Light-adapted fluorescence yield in the dark (`Kp=max`)"
    _f_o′::FT = 0
    "Reversible NPQ rate constant (initially zero)"
    _k_npq_rev::FT = 0
    "Rate constant for photochemistry"
    _k_p::FT = 4
    "Non-Photochemical quenching "
    _npq::FT = 0
    "Energy quenching"
    _q_e::FT = 0
    "Photochemical quenching"
    _q_p::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-18: add the struct of Cytochrome reaction center
#     2022-Feb-07: add more fields
#     2022-Feb-10: add K_X, ϵ_1, and ϵ_2 fields
#     2022-Mar-01: add more fields: η_c and η_l
#     2022-Mar-01: delete Η_C and Η_L, move η_c and η_l to photosynthesis model
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores reaction center information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct CytochromeReactionCenter{FT<:AbstractFloat} <:AbstractReactionCenter{FT}
    # Constant coefficients
    "Fraction of absorbed light used by PSI ETR"
    F_PSI::FT = 0.41 / (0.41 + 0.44)
    "Rate constant of consititutive heat loss from the antennae `[ns⁻¹]`"
    K_D::FT = 0.55
    "Rate constant of fluorescence `[ns⁻¹]`"
    K_F::FT = 0.05
    "Rate constant of photochemistry for PS I `[ns⁻¹]`"
    K_PSI::FT = 14.5
    "Rate constant of photochemistry for PS II `[ns⁻¹]`"
    K_PSII::FT = 4.5
    "Rate constant of excitation sharing for PS II `[ns⁻¹]`"
    K_U::FT = 2
    "Rate constant of regulated heat loss via oxidized PS I center `[s⁻¹]`"
    K_X::FT = 14.5
    "Maximal PS I photochemical yield"
    Φ_PSI_MAX::FT = K_PSI / (K_D + K_F + K_PSI)

    # Diagnostic variables
    "Weight factor that PSI fluorescence reaches sensor (after reabsorption)"
    ϵ_1::FT = 0
    "Weight factor that PSII fluorescence reaches sensor (after reabsorption)"
    ϵ_2::FT = 1
    "Fluorescence yield"
    ϕ_f::FT = 0
    "Photochemical yield"
    ϕ_p::FT = 0
end
