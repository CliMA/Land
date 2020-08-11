###############################################################################
#
# Create RTContainer
#
###############################################################################
"""
    create_rt_container(can::Canopy4RT{FT}, angles::SolarAngles{FT}) where {FT<:AbstractFloat}

Create an [`RTContainer`](@ref), given
- `can` [`Canopy4RT`](@ref) type struct, providing canopy structure information
- `can_opt` [`CanopyOpticals`](@ref) type, where optical properties is stored
- `angles` [`SolarAngles`](@ref) type struct, defining sun-sensor geometry
"""
function create_rt_container(
            can::Canopy4RT{FT},
            can_opt::CanopyOpticals{FT},
            angles::SolarAngles{FT}
) where {FT<:AbstractFloat}
    @unpack lazitab, litab = can;
    @unpack psi = angles;

    cos_ttlo  = cosd.(lazitab);
    cos_philo = cosd.(lazitab .- psi);
    cos_ttli  = cosd.(litab);
    sin_ttli  = sind.(litab);
    vol_scatt = ones(FT, 4);

    _Cs = cos_ttli .* 0; # [nli]
    _Ss = sin_ttli .* 0; # [nli]
    _Co = cos_ttli .* 0; # [nli]
    _So = sin_ttli .* 0; # [nli]
    _1s = ones(FT,1,length(lazitab)); #[1, nlazi]
    cds = _Cs * _1s .+ _Ss * cos_ttlo' ; # [nli, nlazi]
    cdo = _Co * _1s .+ _So * cos_philo'; # [nli, nlazi]
    _2d = _Co * _1s                    ; # [nli, nlazi]

    τ_dd   = can_opt.sigb .* 0;
    τ_sd   = can_opt.sigb .* 0;
    ρ_dd   = can_opt.sigb .* 0;
    ρ_sd   = can_opt.sigb .* 0;
    dnorm  = can_opt.sigb[:,1] .* 0;
    piloc2 = can_opt.sigb .* 0;
    piloc  = can_opt.sigb[:,1] .* 0;
    pilos  = can_opt.sigb[:,1] .* 0;
    pilo   = can_opt.sigb[:,1] .* 0;

    return RTContainer{FT}(cos_ttlo,
                           cos_philo,
                           cos_ttli,
                           sin_ttli,
                           vol_scatt,
                           _Cs,
                           _Ss,
                           _Co,
                           _So,
                           _1s,
                           cds,
                           cdo,
                           _2d,
                           τ_dd,
                           τ_sd,
                           ρ_dd,
                           ρ_sd,
                           dnorm,
                           pilo,
                           piloc2,
                           piloc,
                           pilos)
end
