#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Jan-15: isolate the reaction center from Leaf in Photosynthesis.jl
#
#######################################################################################################################################################################################################
"""
$(TYPEDEF)

Structure that stores reaction center information

# Fields
$(TYPEDFIELDS)
"""
mutable struct PhotosynthesisReactionCenter{FT<:AbstractFloat}
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
    ϕ_s::FT
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
    PhotosynthesisReactionCenter{FT}() where {FT<:AbstractFloat}

Constructor of `PhotosynthesisReactionCenter`
"""
PhotosynthesisReactionCenter{FT}() where {FT<:AbstractFloat} = PhotosynthesisReactionCenter{FT}(0.5, 0.85, 0.05, 4, 4/(0.85+0.05+4), 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0);
