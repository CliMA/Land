###############################################################################
#
# Create SoilOpticals
#
###############################################################################
"""
create_soil_opticals(wl_set::WaveLengths{FT}) where {FT<:AbstractFloat}

Create [`SoilOpticals`](@ref) struct, given
- `wl_set` [`WaveLengths`](@ref) type wave length parameter set
"""
function create_soil_opticals(
        wl_set::WaveLengths{FT}
) where {FT<:AbstractFloat}
    albedo_SW = FT(0.2) * ones(FT, length(wl_set.WL));
    emsvty_SW = 1 .- albedo_SW;

return SoilOpticals{FT}(wl_set.WL,
                        albedo_SW,
                        emsvty_SW,
                        FT[0.1],
                        FT(290.0));;
end
