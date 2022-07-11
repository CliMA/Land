#=
###############################################################################
#
# Diff function to minimize by ConstrainedRootSolvers
# Useful for nocturnal stomatal conductance model
#
###############################################################################
"""
    dRdE(clayer::CanopyLayer{FT},
         envir::AirLayer{FT},
         LAI::FT
    ) where {FT<:AbstractFloat}

Calculate the marginal decrease of respiration rate, given
- `clayer` [`CanopyLayer`](@ref) type of struct
- `envir` `AirLayer` type struct
- `LAI` Total leaf area index of whole canopy
"""
function dRdE(
            clayer::CanopyLayer{FT},
            envir::AirLayer{FT}
) where {FT<:AbstractFloat}
    @unpack T = clayer;

    return dTdE(clayer, envir) * clayer.ps.PSM.r_d * clayer.ps.PSM.TD_R.ΔHA / GAS_R(FT) / T^2
end




"""
    dTdE(clayer::CanopyLayer{FT},
         envir::AirLayer{FT}
    ) where {FT<:AbstractFloat}

Calculate the margian decrease in leaf temperature, given
- `clayer` [`CanopyLayer`](@ref) type of struct
- `envir` `AirLayer` type struct
- `LAI` Total leaf area index of whole canopy
"""
function dTdE(
            clayer::CanopyLayer{FT},
            envir::AirLayer{FT}
) where {FT<:AbstractFloat}
    @unpack T, tLAI, width = clayer;
    @unpack wind = envir;

    lambda = latent_heat_vapor(T) * M_H₂O(FT);
    gbe    = FT(0.189 * sqrt(wind/(0.72*width)));
    emis   = FT(0.97);
    denom  = 2 * CP_D_MOL(FT) * gbe + 4 / tLAI * K_STEFAN(FT) * emis * T^3;
    ∂T∂E   = lambda / denom;

    return ∂T∂E
end




"""
    dΘdE(clayer::CanopyLayer{FT},
         sm::OSMWang{FT}(),
         g_sw::FT
    ) where {FT<:AbstractFloat}

Calculate the margian carbon cost related to nighttime transpiration, given
- `clayer` [`CanopyLayer`](@ref) type of struct
- `sm` [`OSMWang`](@ref) type stomatal model
- `g_sw` Given leaf level stomatal conductance
"""
function dΘdE(
            clayer::CanopyLayer{FT},
            sm::OSMWang{FT},
            g_sw::FT
) where {FT<:AbstractFloat}
    @unpack APAR_m, ec, envir_m, ff = clayer;
    @unpack P_AIR, p_H₂O = envir_m
    clayer.ps_m.apar = APAR_m;

    # calculate g_lc and a_net using memory clayer and envir
    g_sc = g_sw / FT(1.6);
    g_lc = 1 / (1/g_sc + 1/clayer.ps_m.g_CO₂_b);
    g_bw = clayer.ps_m.g_CO₂_b * FT(1.35);
    g_lw = 1 / (1/g_sw + 1/g_bw);
    flow = g_lw * max(1, clayer.ps_m.p_H₂O_sat - p_H₂O) / P_AIR;
    leaf_photosynthesis!(clayer.ps_m, envir_m, GCO₂Mode(), g_lc);

    return clayer.ps_m.PSM.a_net / (ec - flow) * ff
end
=#
