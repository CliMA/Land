###############################################################################
#
# Diff function to minimize by ConstrainedRootSolvers
# Useful for nocturnal stomatal conductance model
#
###############################################################################
"""
    dRdE(photo_set::AbstractPhotoModelParaSet{FT},
         clayer::CanopyLayer{FT},
         envir::AirLayer{FT},
         LAI::FT
    ) where {FT<:AbstractFloat}

Calculate the marginal decrease of respiration rate, given
- `photo_set` `AbstractPhotoModelParaSet` type struct
- `clayer` [`CanopyLayer`](@ref) type of struct
- `envir` `AirLayer` type struct
- `LAI` Total leaf area index of whole canopy
"""
function dRdE(
            photo_set::AbstractPhotoModelParaSet{FT},
            clayer::CanopyLayer{FT},
            envir::AirLayer{FT}
) where {FT<:AbstractFloat}
    @unpack ps, T = clayer;
    @unpack Rd = ps;

    return dTdE(clayer, envir) * Rd * photo_set.ReT.ΔHa_to_RT25 * T_25(FT) /
           T^2
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
    denom  = 2 * CP_D_MOL(FT) * gbe +
             4 / tLAI * K_STEFAN(FT) * emis * T^3;
    ∂T∂E   = lambda / denom;

    return ∂T∂E
end




"""
    dΘdE(photo_set::AbstractPhotoModelParaSet{FT},
         clayer::CanopyLayer{FT},
         sm::OSMWang{FT}(),
         g_sw::FT
    ) where {FT<:AbstractFloat}

Calculate the margian carbon cost related to nighttime transpiration, given
- `photo_set` `AbstractPhotoModelParaSet` type struct
- `clayer` [`CanopyLayer`](@ref) type of struct
- `sm` [`OSMWang`](@ref) type stomatal model
- `g_sw` Given leaf level stomatal conductance
"""
function dΘdE(
            photo_set::AbstractPhotoModelParaSet{FT},
            clayer::CanopyLayer{FT},
            sm::OSMWang{FT},
            g_sw::FT
) where {FT<:AbstractFloat}
    @unpack APAR_m, ec, envir_m, ff, ps_m = clayer;
    @unpack p_atm, p_H₂O = envir_m
    ps_m.APAR = APAR_m;

    # calculate g_lc and a_net using memory clayer and envir
    g_sc = g_sw / FT(1.6);
    g_lc = 1 / (1/g_sc + 1/ps_m.g_bc);
    g_bw = ps_m.g_bc * FT(1.35);
    g_lw = 1 / (1/g_sw + 1/g_bw);
    flow = g_lw * max(1, ps_m.p_sat - p_H₂O) / p_atm;
    leaf_photosynthesis!(photo_set, ps_m, envir_m, GCO₂Mode(), g_lc);

    return ps_m.An / (ec - flow) * ff
end




"""
    nocturnal_diff!(
                x::FT,
                photo_set::AbstractPhotoModelParaSet{FT},
                clayer::CanopyLayer{FT},
                envir::AirLayer{FT},
                sm::OSMWang{FT}()
    ) where {FT<:AbstractFloat}

Calculate the difference between marginal gain and risk of nighttime
    transpiration, given
- `x` Given leaf level stomatal conductance
- `photo_set` `AbstractPhotoModelParaSet` type struct
- `clayer` [`CanopyLayer`](@ref) type of struct
- `envir` `AirLayer` type struct
- `sm` [`OSMWang`](@ref) type stomatal model
"""
function nocturnal_diff!(
            x::FT,
            photo_set::AbstractPhotoModelParaSet{FT},
            clayer::CanopyLayer{FT},
            envir::AirLayer{FT},
            sm::OSMWang{FT}
) where {FT<:AbstractFloat}
    ∂R∂E = dRdE(photo_set, clayer, envir);
    ∂Θ∂E = dΘdE(photo_set, clayer, sm, x);

    return ∂R∂E - ∂Θ∂E
end
