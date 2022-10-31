"""
    spac_beta_max(spac::SPACMono{FT}, beta::BetaGLinearPsoil{FT}) where {FT<:AbstractFloat}

Compute the beta tuning factor for SPAC by taking the maximum, given
- `spac` SPAC
- `beta` `BetaGLinearPsoil` type scheme

"""
function spac_beta_max(spac::SPACMono{FT}, beta::BetaGLinearPsoil{FT}) where {FT<:AbstractFloat}
    _βm::FT = 0;

    for _i in eachindex(spac.plant_hs.roots)
        _β = β_factor(spac.plant_hs.leaves[1], spac.plant_hs.roots[_i].sh, beta, FT(0), spac.plant_hs.roots[_i].p_ups, spac.swc[_i]);
        _βm = max(_β, _βm);
    end;

    return _βm
end
