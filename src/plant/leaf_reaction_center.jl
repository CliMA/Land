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


"""
    VJPReactionCenter{FT}() where {FT<:AbstractFloat}

Constructor of `VJPReactionCenter`
"""
VJPReactionCenter{FT}() where {FT<:AbstractFloat} = VJPReactionCenter{FT}(0.5, 0.85, 0.05, 4, 4/(0.85+0.05+4), 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-18: add the struct of Cytochrome reaction center
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
    "Coupling efficiency of cyclic electron flow `[mol ATP mol⁻¹ e⁻]`"
    Η_C::FT
    "Coupling efficiency of linear electron flow `[mol ATP mol⁻¹ e⁻]`"
    Η_L::FT
    "Maximal PS I photochemical yield"
    Φ_PSI_MAX::FT

    # prognostic variables that change with time

    # dignostic variables that change with time
end
