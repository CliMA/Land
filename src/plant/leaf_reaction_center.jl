#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2022-Jan-14: add van der Tol model struct
#     2022-Jan-24: fix documentation
#     2022-Feb-07: remove the hierachy from abstract fluorescence model
#     2022-Feb-07: move struct definition as a field of VJPReactionCenter
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
struct VanDerTolFluorescenceModel{FT<:AbstractFloat}
    "Fitting parameter K_0"
    K_0::FT
    "Fitting parameter α"
    K_A::FT
    "Fitting parameter β"
    K_B::FT
end


#######################################################################################################################################################################################################
#
# Changes to the constructor
# General
#     2022-Jan-14: migrate from Photosynthesis.jl
#     2022-Feb-24: set default to true to make the response more reasonable
# To do
#     TODO: examine why van der Tol et al has the nonstressed parameter set that are so far away from the stressed one
# Sources
#     van der Tol et al. (2013) Models of fluorescence and photosynthesis for interpreting measurements of solar-induced chlorophyll fluorescence
#
#######################################################################################################################################################################################################
"""

    VanDerTolFluorescenceModel{FT}(drought::Bool = true) where {FT<:AbstractFloat}

Constructor for `VanDerTolFluorescenceModel` fluorescence model, given
- `drought` If true, use parameters trained from drought stressed plant. Default is `true`.

---
# Examples
```julia
vdt = VanDerTolFluorescenceModel{Float64}();
vdt = VanDerTolFluorescenceModel{Float64}(false);
```
"""
VanDerTolFluorescenceModel{FT}(drought::Bool = true) where {FT<:AbstractFloat} = (
    if drought
        return VanDerTolFluorescenceModel{FT}(5.01, 1.93, 10)
    else
        return VanDerTolFluorescenceModel{FT}(2.48, 2.83, 0.114)
    end
);


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jan-24: abstractize the reaction center
#     2022-Jan-25: fix documentation
#     2022-Feb-07: fix documentation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores reaction center information

Hierachy of the `AbstractReactionCenter`
- [`VJPReactionCenter`](@ref)
- [`CytochromeReactionCenter`](@ref)
"""
abstract type AbstractReactionCenter{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-15: isolate the reaction center from Leaf in Photosynthesis.jl
#     2022-Jan-25: fix documentation
#     2022-Feb-07: add fluorescence model as a field
#     2022-Mar-01: fix documentation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores reaction center information

# Fields

$(TYPEDFIELDS)

"""
mutable struct VJPReactionCenter{FT<:AbstractFloat} <:AbstractReactionCenter{FT}
    # parameters that do not change with time
    "Fraction of absorbed light used by PSII ETR"
    F_PSII::FT
    "Fluorescence model"
    FLM::VanDerTolFluorescenceModel{FT}
    "Rate constant for thermal dissipation"
    K_D::FT
    "Rate constant for fluorescence"
    K_F::FT
    "Maximal rate constant for photochemistry"
    K_P_MAX::FT
    "max PSII yield (k_npq_r=0, all RC open)"
    Φ_PSII_MAX::FT

    # prognostic variables that change with time
    "Reversible NPQ rate constant (initially zero)"
    k_npq_rev::FT
    "Sustained NPQ rate constant (for seasonal changes, default is zero)"
    k_npq_sus::FT
    "Rate constant for photochemistry"
    k_p::FT
    "Non-Photochemical quenching "
    npq::FT
    "Fluorescence yield"
    ϕ_f::FT
    "Photochemical yield"
    ϕ_p::FT

    # dignostic variables that change with time
    "Dark adapted yield (`Kp=0`)"
    f_m::FT
    "Light adapted yield (`Kp=0`)"
    f_m′::FT
    "Dark-adapted fluorescence yield (`Kp=max`)"
    f_o::FT
    "Light-adapted fluorescence yield in the dark (`Kp=max`)"
    f_o′::FT
    "Energy quenching"
    q_e::FT
    "Photochemical quenching"
    q_p::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jan-15: isolate the reaction center from Leaf in Photosynthesis.jl
#     2022-Jan-25: fix documentation
#
#######################################################################################################################################################################################################
"""

    VJPReactionCenter{FT}() where {FT<:AbstractFloat}

Constructor of `VJPReactionCenter`

---
# Examples
```julia
rc = VJPReactionCenter{Float64}();
```
"""
VJPReactionCenter{FT}() where {FT<:AbstractFloat} = (
    return VJPReactionCenter{FT}(
                0.5,                                # F_PSII
                VanDerTolFluorescenceModel{FT}(),   # FLM
                0.85,                               # K_D
                0.05,                               # K_F
                4,                                  # K_P_MAX
                4/(0.85+0.05+4),                    # Φ_PSII_MAX
                0,                                  # k_npq_rev
                0,                                  # k_npq_sus
                4,                                  # k_p
                0,                                  # npq
                0,                                  # ϕ_f
                0,                                  # ϕ_p
                0,                                  # f_m
                0,                                  # f_m′
                0,                                  # f_o
                0,                                  # f_o′
                0,                                  # q_e
                0)                                  # q_p
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-18: add the struct of Cytochrome reaction center
#     2022-Jan-25: fix documentation
#     2022-Feb-07: add more fields to work with Photosynthesis v0.3.1
#     2022-Feb-10: add K_X, ϵ_1, and ϵ_2 fields
#     2022-Mar-01: fix documentation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores reaction center information

# Fields

$(TYPEDFIELDS)

"""
mutable struct CytochromeReactionCenter{FT<:AbstractFloat} <:AbstractReactionCenter{FT}
    # parameters that do not change with time
    "Fraction of absorbed light used by PSI ETR"
    F_PSI::FT
    "Rate constant of consititutive heat loss from the antennae `[ns⁻¹]`"
    K_D::FT
    "Rate constant of fluorescence `[ns⁻¹]`"
    K_F::FT
    "Rate constant of photochemistry for PS I `[ns⁻¹]`"
    K_PSI::FT
    "Rate constant of photochemistry for PS II `[ns⁻¹]`"
    K_PSII::FT
    "Rate constant of excitation sharing for PS II `[ns⁻¹]`"
    K_U::FT
    "Rate constant of regulated heat loss via oxidized PS I center `[s⁻¹]`"
    K_X::FT
    "Coupling efficiency of cyclic electron flow `[mol ATP mol⁻¹ e⁻]`"
    Η_C::FT
    "Coupling efficiency of linear electron flow `[mol ATP mol⁻¹ e⁻]`"
    Η_L::FT
    "Maximal PS I photochemical yield"
    Φ_PSI_MAX::FT

    # prognostic variables that change with time
    "Weight factor that PSI fluorescence reaches sensor (after reabsorption)"
    ϵ_1::FT
    "Weight factor that PSII fluorescence reaches sensor (after reabsorption)"
    ϵ_2::FT
    "Coupling efficiency of cyclic electron flow `[mol ATP mol⁻¹ e⁻]`"
    η_c::FT
    "Coupling efficiency of linear electron flow `[mol ATP mol⁻¹ e⁻]`"
    η_l::FT
    "Fluorescence yield"
    ϕ_f::FT
    "Photochemical yield"
    ϕ_p::FT

    # dignostic variables that change with time
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jan-24: add constructor for CytochromeReactionCenter
#     2022-Jan-25: fix documentation
#     2022-Feb-07: make the constructor more readable
#     2022-Feb-10: add K_X, ϵ_1, and ϵ_2 fields
#
#######################################################################################################################################################################################################
"""
    CytochromeReactionCenter{FT}() where {FT<:AbstractFloat}

Constructor of `CytochromeReactionCenter`

---
# Examples
```julia
rc = CytochromeReactionCenter{Float64}();
```
"""
CytochromeReactionCenter{FT}() where {FT<:AbstractFloat} = (
    return CytochromeReactionCenter{FT}(
                0.41/(0.41+0.44),           # F_PSI
                0.55,                       # K_D
                0.05,                       # K_F
                14.5,                       # K_PSI
                4.5,                        # K_PSII
                0,                          # K_U
                14.5,                       # K_X
                1,                          # Η_C
                0.75,                       # Η_L
                14.5 / (14.5+0.55+0.05),    # Φ_PSI_MAX
                0,                          # ϵ_1
                1,                          # ϵ_2
                1,                          # η_c
                0.75,                       # η_l
                0,                          # ϕ_f
                0)                          # ϕ_p
);
