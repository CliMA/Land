###############################################################################
#
# Update leaf sensible and latent heat flux
#
###############################################################################
"""
    leaf_heat_flux!(leaf::Leaf{FT}, envir::AirLayer{FT})

Update leaf sensible and latent fluexes, given
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct

The sensible heat flux is computed using
```math
H = 2 ⋅ ρ ⋅ c_p ⋅ g_b ⋅ (T - T_{air})
```
where 2 describes the two-sided leaf, ρ is the density of dry air
(``ρ = \\dfrac{p_{atm} ⋅ M_{dryair}}{R ⋅ T_{air}}``, ``M_{dryair}`` is the
molar mass of dry air and ``R`` is the gas constant), ``c_p`` is heat capacity
of dry air at constant pressure, ``g_b`` is the leaf boundary layer conductance
for heat transfer in `[m s⁻¹]` (`[mol m⁻² s⁻¹]` times
``\\dfrac{RT_{air}}{P_{atm}}`` gives `[m s⁻¹]` ), ``T`` is leaf temperature,
and ``T_{air}`` is air temperature. Thus, the formulation is simplified to
```math
H = 2 ⋅ c_p ⋅ g_{bw} ⋅ (T - T_{air}) ⋅ M_{dryair}
```

The latent heat flux is computed using
```math
LE = LV ⋅ E
```
where ``LV`` is specific latent heat of evporation in `[J mol⁻¹]`, and ``E`` is
evaporation rate in `[mol m⁻² s⁻¹]`.
"""
function leaf_heat_flux!(
            leaf::Leaf{FT},
            envir::AirLayer{FT}
            ) where {FT<:AbstractFloat}
    @unpack g_bw, T = leaf;
    @unpack t_air = envir;

    leaf.H  = 2 * FT(CP_D) * g_bw * (T - t_air) * FT(MOLMASS_DRYAIR);
    leaf.LE = leaf.LV * leaf.e;

    return nothing
end

function leaf_heat_flux!(
            leaf::Leaves{FT},
            envir::AirLayer{FT}
            ) where {FT<:AbstractFloat}
    @unpack g_bw, T = leaf;
    @unpack t_air = envir;

    leaf.H  .= 2 .* FT(CP_D) .* g_bw .* (T - t_air) .* FT(MOLMASS_DRYAIR);
    leaf.LE .= leaf.LV .* leaf.e;

    return nothing
end
