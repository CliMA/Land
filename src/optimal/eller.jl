###############################################################################
#
# Diff function to minimize by ConstrainedRootSolvers
#
###############################################################################
"""
    envir_diff!(x::FT, photo_set::AbstractPhotoModelParaSet{FT}, canopyi::CanopyLayer{FT}, hs::LeafHydraulics{FT}, envir::AirLayer{FT}, sm::AbstractStomatalModel{FT}, ind::Int) where {FT<:AbstractFloat}

Calculate the difference to be minimized for a given
- `x` Assumed leaf diffusive conductance
- `photo_set`[`C3ParaSet`] or [`C4ParaSet`] type parameter set
- `canopyi`[`CanopyLayer`](@ref) type struct
- `hs` Leaf hydraulic system
- `envir`[`AirLayer`] type struct
- `sm` Stomatal model option (photo_set.Sto)
- `int` Nth leaf in the canopy layer
"""
function envir_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet,
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            hs::LeafHydraulics{FT},
            sm::ESMGentine{FT},
            ind::Int
            ) where {FT<:AbstractFloat}
    # unpack variables
    @unpack p_sat, ps = canopyi;
    @unpack p_atm, p_H₂O = envir;
    g_bc  = canopyi.g_bc[ind];
    g_bw  = canopyi.g_bw[ind];
    g_m   = canopyi.g_m[ind];

    # update photosynthesis for ps
    ps.g_lc = x;
    leaf_photo_from_glc!(photo_set, ps, envir);

    # calculate flow rate and β
    g_sc = 1 / max( 1/x - 1/g_m - 1/g_bc, FT(1.6e-3) );
    g_sw = g_sc * FT(1.6);
    g_lw = 1 / (1/g_sw + 1/g_bw);
    e_lf = g_lw * (p_sat - p_H₂O) / p_atm;
    k_lf = leaf_xylem_risk(hs, e_lf);

    # calculate g_sw from stomatal model
    g_md = empirical_gsw_from_model(sm, ps, envir, k_lf);
    g_md = min(canopyi.g_max, g_md);

    # calculate model predicted g_lc
    g_lm = 1 / (FT(1.6)/g_md + 1/g_bc + 1/g_m);

    return g_lm - x
end