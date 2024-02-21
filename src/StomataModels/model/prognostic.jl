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
                sm::OSMWang{FT},
                Δt::FT
    ) where {FT<:AbstractFloat}

Update g_sw prognostically, given
- `clayer` A `CanopyLayer` type struct
- `envir` `AirLayer` type environmental conditions
- `sm` `EmpiricalStomatalModel` or `OSMWang` type stomatal model
- `β` Tune factor to stomatal g1. 1 for AbstractBetaV mode
- `Δt` Time interval for prognostic stomatal conductance
- `photo_set` `AbstractPhotoModelParaSet` type photosynthesis model, currently
    supports [`OSMWang`](@ref) model only
- `hs` Leaf hydraulic system
"""
function prognostic_gsw! end

prognostic_gsw!(clayer::CanopyLayer{FT}, envir::AirLayer{FT}, sm::EmpiricalStomatalModel{FT}, β::FT, Δt::FT) where {FT<:AbstractFloat} = (
    # unpack values
    (; g_bc, g_bw, g_lc, g_lw, g_m, g_sc, g_sw, n_leaf) = clayer;

    # update g_sw
    for iLF in 1:n_leaf
        gsw_ss = max(0, stomatal_conductance(sm, clayer, envir, β, iLF));
        g_sw[iLF] += (gsw_ss - g_sw[iLF]) / clayer.τ_esm * Δt;

        # update g_lw, gsc, and g_lc as well
        g_lw[iLF] = 1 / ( 1/g_sw[iLF]  + 1/g_bw[iLF] );
        g_sc[iLF] = g_sw[iLF] / FT(1.6);
        g_lc[iLF] = 1 / ( 1/g_sc[iLF] + 1/g_m[iLF] + 1/g_bc[iLF] );
    end;

    return nothing
);

prognostic_gsw!(photo_set::AbstractPhotoModelParaSet{FT}, clayer::CanopyLayer{FT}, hs::LeafHydraulics{FT}, envir::AirLayer{FT}, sm::OSMWang{FT}, Δt::FT) where {FT<:AbstractFloat} = (
    # unpack values
    (; APAR, g_bc, g_bw, g_lc, g_lw, g_m, g_sc, g_sw, n_leaf) = clayer;

    # update g_sw
    for iLF in 1:n_leaf
        clayer.ps.APAR = clayer.APAR[iLF];
        leaf_ETR!(photo_set, clayer.ps);

        if APAR[iLF] > 1
            diff = solution_diff!(g_sw[iLF], photo_set, clayer, hs, envir, sm, GswDrive(), iLF);
            g_sw[iLF] += diff * clayer.τ_osm * Δt;
        else
            diff = nocturnal_diff!(g_sw[iLF], photo_set, clayer, envir, sm);
            g_sw[iLF] += diff * clayer.τ_noc * Δt;
        end

        # update g_lw, gsc, and g_lc as well
        g_lw[iLF] = 1 / ( 1/g_sw[iLF]  + 1/g_bw[iLF] );
        g_sc[iLF] = g_sw[iLF] / FT(1.6);
        g_lc[iLF] = 1 / ( 1/g_sc[iLF] + 1/g_m[iLF] + 1/g_bc[iLF] );
    end;

    return nothing
);
