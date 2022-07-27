###############################################################################
#
# Calculate zenith angle
#
###############################################################################
"""
    zenith_angle(
                latd::FT,
                decd::FT,
                lhad::FT
    ) where {FT<:AbstractFloat}
    zenith_angle(
                latd::FT,
                day::FT,
                hour::FT,
                minute::FT
    ) where {FT<:AbstractFloat}

Calculate the zenith angle, given
- `latd` Latitude in degree
- `decd` Declination of the Sun in degree
- `lhad` Local hour angle in degree
- `day` Day of year
- `hour` Hour of day
- `minute` Minute of hour
"""
function zenith_angle(
            latd::FT,
            decd::FT,
            lhad::FT
) where {FT<:AbstractFloat}
    cosz = sind(latd) * sind(decd) + cosd(latd) * cosd(decd) * cosd(lhad);

    return acosd(cosz)
end




function zenith_angle(
            latd::FT,
            day::FT,
            hour::FT,
            minute::FT
) where {FT<:AbstractFloat}
    _deg::FT = 360 / YEAR_D(FT) * (day + (hour+minute/60) / 24 + 10);
    decd::FT = -FT(23.44) * cosd(_deg);
    lhad::FT = (hour-12) * 15;

    return zenith_angle(latd, decd, lhad)
end




function zenith_angle(
            latd::FT,
            toy::FT
) where {FT<:AbstractFloat}
    _deg::FT = 360 / YEAR_D(FT) * (toy + 10);
    decd::FT = -FT(23.44) * cosd(_deg);
    lhad::FT = ((toy%1) - FT(0.5)) * 360;

    return zenith_angle(latd, decd, lhad)
end
