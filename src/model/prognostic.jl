###############################################################################
#
# Update g_sw using different stomatal models
#
###############################################################################
"""
    prognostic_gsw!(
                clayer::CanopyLayer{FT},
                envir::AirLayer{FT},
                sm::EmpiricalStomatalModel{FT},
                β::FT,
                Δt::FT
    ) where {FT<:AbstractFloat}
    prognostic_gsw!(
                photo_set::AbstractPhotoModelParaSet{FT},
                clayer::CanopyLayer{FT},
                hs::LeafHydraulics{FT},
                envir::AirLayer{FT},
                sm::OptimizationStomatalModel{FT},
                Δt::FT
    ) where {FT<:AbstractFloat}

Update g_sw prognostically, given
- `clayer` A `CanopyLayer` type struct
- `envir` `AirLayer` type environmental conditions
- `sm` `EmpiricalStomatalModel` or `OSMWang` type stomatal model
- `β` Tune factor to stomatal g1. 1 for AbstractBetaV mode
- `Δt` Time interval for prognostic stomatal conductance
- `photo_set` `AbstractPhotoModelParaSet` type photosynthesis model
- `hs` Leaf hydraulic system
"""
function prognostic_gsw!(
            clayer::CanopyLayer{FT},
            envir::AirLayer{FT},
            sm::EmpiricalStomatalModel{FT},
            β::FT,
            Δt::FT
) where {FT<:AbstractFloat}
    # unpack values
    @unpack g_sw, n_leaf = clayer;

    # update g_sw
    for iLF in 1:n_leaf
        gsw_ss = max(sm.g0,
                     stomatal_conductance(sm, clayer, envir, β, iLF));
        g_sw[iLF] += (gsw_ss - g_sw[iLF]) / clayer.τ_esm * Δt;
    end

    return nothing
end




function prognostic_gsw!(
            photo_set::AbstractPhotoModelParaSet{FT},
            clayer::CanopyLayer{FT},
            hs::LeafHydraulics{FT},
            envir::AirLayer{FT},
            sm::OptimizationStomatalModel{FT},
            Δt::FT
) where {FT<:AbstractFloat}
    # unpack values
    @unpack An, ec, g_bc, g_bw, g_lw, g_m, g_sw, n_leaf, p_sat, ps = clayer;
    @unpack p_atm, p_H₂O = envir;

    # update g_sw
    for iLF in 1:n_leaf
        diff = solution_diff!(g_sw[iLF], photo_set, clayer, hs, envir, sm,
                              GswDrive(), iLF);
        g_sw[iLF] += diff * clayer.τ_osm * Δt;
    end

    return nothing
end
