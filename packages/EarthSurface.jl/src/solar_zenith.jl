#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Sep-08: move function from SoilPlantAirContinuum.jl
#     2022-Sep-08: make YEAR_D an option so that one can use ClimaCache value to replace it
#     2022-Oct-19: use YEAR_D function from EmeraldConstants
#
#######################################################################################################################################################################################################
"""

    solar_zenith_angle(lat::FT, decl::FT, lha::FT) where {FT<:AbstractFloat}
    solar_zenith_angle(lat::FT, day::FT, hour::FT, minute::FT) where {FT<:AbstractFloat}
    solar_zenith_angle(lat::FT, fdoy::FT) where {FT<:AbstractFloat}

Return the solar zenith angle, given
- `lat` Latitude in °
- `decl` Declination of the Sun in °
- `lha` Local hour angle in °
- `day` Day of year
- `hour` Hour of day
- `minute` Minute of hour
- `fdoy` Day of year (digits after decimal for time of day)

"""
function solar_zenith_angle end

solar_zenith_angle(lat::FT, decl::FT, lha::FT) where {FT<:AbstractFloat} = (
    _cos_sza = sind(lat) * sind(decl) + cosd(lat) * cosd(decl) * cosd(lha);

    return acosd(_cos_sza)
);

solar_zenith_angle(lat::FT, day::FT, hour::FT, minute::FT) where {FT<:AbstractFloat} = (
    _deg  = 360 / YEAR_D(FT) * (day + (hour + minute / 60) / 24 + 10);
    _decl = -FT(23.44) * cosd(_deg);
    _lha  = (hour - 12) * 15;

    return solar_zenith_angle(lat, _decl, _lha)
);

solar_zenith_angle(lat::FT, fdoy::FT) where {FT<:AbstractFloat} = (
    _deg  = 360 / YEAR_D(FT) * (fdoy + 10);
    _decl = -FT(23.44) * cosd(_deg);
    _lha  = ((fdoy % 1) - FT(0.5)) * 360;

    return solar_zenith_angle(lat, _decl, _lha)
);
