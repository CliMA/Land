#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jan-24: abstractize the reaction center
#     2022-Jan-25: fix documentation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores reaction center information

Hierachy of the `AbstractFluorescenceModel`
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
    "Rate constant for photochemistry (all reaction centers open)"
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
VJPReactionCenter{FT}() where {FT<:AbstractFloat} = VJPReactionCenter{FT}(0.5, 0.85, 0.05, 4, 4/(0.85+0.05+4), 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-18: add the struct of Cytochrome reaction center
#     2022-Jan-25: fix documentation
#     2022-Feb-07: add more fields to work with Photosynthesis v0.3.1
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
    "Rate constant of consititutive heat loss from the antennae `[s⁻¹]`"
    K_D::FT
    "Rate constant of fluorescence `[s⁻¹]`"
    K_F::FT
    "Rate constant of photochemistry for PS I `[s⁻¹]`"
    K_PSI::FT
    "Rate constant of photochemistry for PS II `[s⁻¹]`"
    K_PSII::FT
    "Rate constant of excitation sharing for PS II `[s⁻¹]`"
    K_U::FT
    "Coupling efficiency of cyclic electron flow `[mol ATP mol⁻¹ e⁻]`"
    Η_C::FT
    "Coupling efficiency of linear electron flow `[mol ATP mol⁻¹ e⁻]`"
    Η_L::FT
    "Maximal PS I photochemical yield"
    Φ_PSI_MAX::FT

    # prognostic variables that change with time
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
                0.41/(0.41+0.44),               # F_PSI
                5.5e8,                          # K_D
                5e7,                            # K_F
                1.45e10,                        # K_PSI
                4.5e9,                          # K_PSII
                0,                              # K_U
                1,                              # Η_C
                0.75,                           # Η_L
                1.45e10 / (1.45e10+5.5e8+5e7),  # Φ_PSI_MAX
                0,                              # ϕ_f
                0)                              # ϕ_p
);
