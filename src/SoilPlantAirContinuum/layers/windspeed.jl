###############################################################################
#
# Surface roughness
#
###############################################################################
"""
    aerodynamic_roughness(h_can::FT) where {FT<:AbstractFloat}

Calculate the aerodynamic roughness for wind, given
- `h_can` Average canopy height
"""
function aerodynamic_roughness(h_can::FT) where {FT<:AbstractFloat}
    return h_can * FT(0.07)
end








###############################################################################
#
# Zero plane displacement
#
###############################################################################
"""
    zero_plane_displacement(h_can::FT) where {FT<:AbstractFloat}

Calculate the zero plane displacement height, given
- `h_can` Average canopy height
"""
function zero_plane_displacement(h_can::FT) where {FT<:AbstractFloat}
    return h_can * 2/3
end








###############################################################################
#
# Calculate wind speed at different height
#
###############################################################################
"""
    wind_speed(us::FT, z0::FT, d::FT, z::FT) where {FT<:AbstractFloat}
    wind_speed(wr::FT, zr::FT, z0::FT, d::FT, z::FT) where {FT<:AbstractFloat}

Calculate the wind speed at height `z`, given
- `us` Friction velocity `[m s⁻¹]`
- `d` Zero plane displacement `[m]`
- `z0` Aerodynamic roughnes `[m]`
- `z` Height `[m]`
- `wr` Wind speed at reference height `[m s⁻¹]`
- `zr` Reference height `[m s⁻¹]`

The equation used is
```math
u(z) = \\dfrac{us}{κ} ⋅ \\log \\left( \\dfrac{z - d}{z0} \\right)
```
where `κ` is the Von Kármán constant.
"""
function wind_speed(
            us::FT,
            z0::FT,
            d::FT,
            z::FT
) where {FT<:AbstractFloat}
    return us / K_VON_KARMAN(FT) * log((z-d)/z0)
end




function wind_speed(
            wr::FT,
            zr::FT,
            z0::FT,
            d::FT,
            z::FT
) where {FT<:AbstractFloat}
    return wr * log((z-d)/z0) / log((zr-d)/z0)
end








###############################################################################
#
# Update wind speed at different canopy layer
#
###############################################################################
"""
    wind_speed!(spac::SPACMono{FT}, us::FT) where {FT<:AbstractFloat}
    wind_speed!(spac::SPACMono{FT}, wr::FT, zr::FT) where {FT<:AbstractFloat}

Update wind speed profile in the canopy, given
- `spac` [`SPACMono`](@ref) type struct
- `us` Friction velocity `[m s⁻¹]`
- `wr` Wind speed at reference height `[m s⁻¹]`
- `zr` Reference height `[m s⁻¹]`
"""
function wind_speed!(
            spac::SPACMono{FT},
            us::FT
) where {FT<:AbstractFloat}
    @unpack air_bounds, n_canopy, wind_d, wind_z0, wind_zs = spac;
    spac.winds .= wind_speed.(us, wind_z0, wind_d, wind_zs);

    return nothing
end




function wind_speed!(
            spac::SPACMono{FT},
            wr::FT,
            zr::FT
) where {FT<:AbstractFloat}
    @unpack air_bounds, n_canopy, wind_d, wind_z0, wind_zs = spac;
    spac.winds .= wind_speed.(wr, zr, wind_z0, wind_d, wind_zs);

    return nothing
end
