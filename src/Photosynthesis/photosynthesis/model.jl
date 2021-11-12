###############################################################################
#
# Calculate photosynthesis from CO₂ partial pressure
#
###############################################################################
"""
    leaf_photosynthesis!(
                photo_set::AbstractPhotoModelParaSet{FT},
                leaf::Leaf{FT},
                envir::AirLayer{FT},
                mode::PCO₂Mode
    ) where {FT<:AbstractFloat}
    leaf_photosynthesis!(
                photo_set::AbstractPhotoModelParaSet{FT},
                leaf::Leaf{FT},
                envir::AirLayer{FT},
                mode::PCO₂Mode,
                p_i::FT
    ) where {FT<:AbstractFloat}
    leaf_photosynthesis!(
                photo_set::C3ParaSet{FT},
                leaf::Leaf{FT},
                envir::AirLayer{FT},
                mode::GCO₂Mode
    ) where {FT<:AbstractFloat}
    leaf_photosynthesis!(
                photo_set::C4ParaSet{FT},
                leaf::Leaf{FT},
                envir::AirLayer{FT},
                mode::GCO₂Mode
    ) where {FT<:AbstractFloat}
    leaf_photosynthesis!(
                photo_set::AbstractPhotoModelParaSet{FT},
                leaf::Leaf{FT},
                envir::AirLayer{FT},
                mode::GCO₂Mode,
                g_lc::FT
    ) where {FT<:AbstractFloat}

Compute leaf photosynthetic rates, given
- `photo_set` [`AbstractPhotoModelParaSet`](@ref) type parameter set
- `leaf` [`Leaf`](@ref) type struct
- `p_i` Given leaf internal CO₂
- `envir` [`AirLayer`](@ref) type struct
- `g_lc` Given leaf diffusive conductance to CO₂

The C3 photosynthesis model is from Farquhar et al. (1980) "A biochemical model
    of photosynthetic CO₂ assimilation in leaves of C3 species."

The C4 photosynthesis model is adapted from Collatz et al. (1992) "Coupled
    photosynthesis-stomatal conductance model for leaves of C4 plants."
"""
function leaf_photosynthesis!(
            photo_set::C3Cytochrome{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            mode::PCO₂Mode
) where {FT<:AbstractFloat}
    leaf_temperature_dependence!(photo_set, leaf, envir);
    leaf_ETR!(photo_set, leaf);
    light_limited_rate!(photo_set, leaf);
    rubisco_limited_rate!(photo_set, leaf);
    product_limited_rate!(photo_set, leaf);
    leaf.Ag = min(leaf.Ac, leaf.Aj, leaf.Ap);
    leaf.An = leaf.Ag - leaf.Rd;
    leaf.J_P680_a = min(leaf.J_P680_c, leaf.J_P680_j, leaf.J_P680_p);
    leaf.J_P700_a = min(leaf.J_P700_c, leaf.J_P700_j, leaf.J_P700_p);

    return nothing
end




function leaf_photosynthesis!(
            photo_set::AbstractPhotoModelParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            mode::PCO₂Mode
) where {FT<:AbstractFloat}
    leaf_temperature_dependence!(photo_set, leaf, envir);
    leaf_ETR!(photo_set, leaf);
    light_limited_rate!(photo_set, leaf);
    rubisco_limited_rate!(photo_set, leaf);
    product_limited_rate!(photo_set, leaf);
    leaf.Ag = min(leaf.Ac, leaf.Aj, leaf.Ap);
    leaf.An = leaf.Ag - leaf.Rd;

    return nothing
end




function leaf_photosynthesis!(
            photo_set::AbstractPhotoModelParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            mode::PCO₂Mode,
            p_i::FT
) where {FT<:AbstractFloat}
    leaf.p_i = p_i;
    leaf_photosynthesis!(photo_set, leaf, envir, mode);

    return nothing
end




function leaf_photosynthesis!(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            mode::GCO₂Mode
) where {FT<:AbstractFloat}
    leaf_temperature_dependence!(photo_set, leaf, envir);
    leaf_ETR!(photo_set, leaf);
    light_limited_rate!(photo_set, leaf, envir);
    rubisco_limited_rate!(photo_set, leaf, envir);
    product_limited_rate!(photo_set, leaf);
    leaf.Ag  = min(leaf.Ac, leaf.Aj, leaf.Ap);
    leaf.An  = leaf.Ag - leaf.Rd;
    leaf.p_i = envir.p_a - leaf.An / leaf.g_lc * envir.p_atm * FT(1e-6);
    leaf.p_s = envir.p_a - leaf.An / leaf.g_bc * envir.p_atm * FT(1e-6);

    return nothing
end




function leaf_photosynthesis!(
            photo_set::C4ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            mode::GCO₂Mode
) where {FT<:AbstractFloat}
    leaf_temperature_dependence!(photo_set, leaf, envir);
    leaf_ETR!(photo_set, leaf);
    light_limited_rate!(photo_set, leaf);
    rubisco_limited_rate!(photo_set, leaf);
    product_limited_rate!(photo_set, leaf, envir);
    leaf.Ag  = min(leaf.Ac, leaf.Aj, leaf.Ap);
    leaf.An  = leaf.Ag - leaf.Rd;
    leaf.p_i = envir.p_a - leaf.An / leaf.g_lc * envir.p_atm * FT(1e-6);
    leaf.p_s = envir.p_a - leaf.An / leaf.g_bc * envir.p_atm * FT(1e-6);

    return nothing
end




function leaf_photosynthesis!(
            photo_set::AbstractPhotoModelParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            mode::GCO₂Mode,
            g_lc::FT
) where {FT<:AbstractFloat}
    leaf.g_lc = g_lc;
    leaf_photosynthesis!(photo_set, leaf, envir, mode);

    return nothing
end
