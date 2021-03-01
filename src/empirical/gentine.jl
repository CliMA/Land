###############################################################################
#
# Gentine model
#
###############################################################################
function stomatal_conductance(
            model::ESMGentine{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            β::FT
) where {FT<:AbstractFloat}
    @unpack g0, g1  = model;
    @unpack An, p_i = leaf;
    @unpack p_atm   = envir;

    return g0 + g1 * p_atm * FT(1e-6) * β * An / p_i
end




function stomatal_conductance(
            model::ESMGentine{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            β::FT
) where {FT<:AbstractFloat}
    @unpack g0, g1  = model;
    @unpack An, p_i = canopyi;
    @unpack p_atm   = envir;

    return g0 .+ g1 * p_atm * FT(1e-6) * β .* An ./ p_i
end




function stomatal_conductance(
            model::ESMGentine{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            β::FT,
            ind::Int
) where {FT<:AbstractFloat}
    @unpack g0, g1  = model;
    @unpack An, p_i = canopyi;
    @unpack p_atm   = envir;

    return g0 + g1 * p_atm * FT(1e-6) * β * An[ind] / p_i[ind]
end








###############################################################################
#
# Diff function to minimize by ConstrainedRootSolvers
# Description in general section
#
###############################################################################
function envir_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet{FT},
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
    k_lf = xylem_risk(hs, e_lf);

    # calculate g_sw from stomatal model
    g_md = stomatal_conductance(sm, ps, envir, k_lf);
    g_md = min(canopyi.g_max, g_md);

    # calculate model predicted g_lc
    g_lm = 1 / (FT(1.6)/g_md + 1/g_bc + 1/g_m);

    return g_lm - x
end
