###############################################################################
#
# Initialize SoilOpticals
#
###############################################################################
"""
    create_soil_opticals(wls::WaveLengths{FT}) where {FT<:AbstractFloat}

Create [`SoilOpticals`](@ref) struct, given
- `wls` [`WaveLengths`](@ref) type struct
"""
function create_soil_opticals(wls::WaveLengths{FT}) where {FT<:AbstractFloat}
    @unpack iWLF, nWL = wls;

    albedo_SW     = FT(0.2) * ones(FT, nWL);
    albedo_SW_SIF = albedo_SW[iWLF];
    emsvty_SW     = 1 .- albedo_SW;

return SoilOpticals{FT}(albedo_SW,
                        albedo_SW_SIF,
                        emsvty_SW,
                        FT[0.1],
                        FT(290.0))
end
