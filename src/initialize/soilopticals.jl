###############################################################################
#
# Create SoilOpticals
#
###############################################################################
"""
    create_soil_opticals(wl_set::WaveLengths{FT}) where {FT<:AbstractFloat}

Create [`SoilOpticals`](@ref) struct, given
- `rt_dim` [`WaveLengths`](@ref) type wave length parameter set
"""
function create_soil_opticals(
            FT,
            rt_dim::RTDimentions
)
    albedo_SW = FT(0.2) * ones(FT, rt_dim.nWL);
    emsvty_SW = 1 .- albedo_SW;

return SoilOpticals{FT}(albedo_SW,
                        emsvty_SW,
                        FT[0.1],
                        FT(290.0));;
end
