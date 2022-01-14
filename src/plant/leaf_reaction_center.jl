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
    k_npq_r::FT
    "Sustained NPQ rate constant (for seasonal changes, default is zero)"
    k_npq_s::FT
    "Rate constant for photochemistry (all reaction centers open)"
    k_p::FT
end


"""
    PhotosynthesisReactionCenter{FT}() where {FT<:AbstractFloat}

Constructor of `PhotosynthesisReactionCenter`
"""
PhotosynthesisReactionCenter{FT}() where {FT<:AbstractFloat} = PhotosynthesisReactionCenter{FT}(0.5, 0.85, 0.05, 4, 4/(0.85+0.05+4), 0, 0, 4);
