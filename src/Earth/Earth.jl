#module Earth




# export public functions
export atmospheric_pressure,
       atmospheric_pressure_ratio,
       ppm_to_Pa,
       zenith_angle




###############################################################################
#
# Atmospheric pressure and partial pressure
#
###############################################################################
const GM2RT25 = GRAVITY * 0.02896968 / RK_25;

"""
    atmospheric_pressure(h::FT)

Calculate the atmospheric pressure, given
- `h` elevation in `[m]`
"""
function atmospheric_pressure_ratio(
            h::FT
            ) where {FT<:AbstractFloat}
    return exp( -h * FT(GM2RT25) );
end




"""
    atmospheric_pressure(h::FT)

Calculate the atmospheric pressure, given
- `h` elevation in `[m]`
"""
function atmospheric_pressure(
            h::FT
            ) where {FT<:AbstractFloat}
    return FT(P_ATM) * atmospheric_pressure_ratio(h);
end




"""
    ppm_to_Pa(h::FT)

Convert ppm to Pa, given
- `h` elevation in `[m]`
"""
function ppm_to_Pa(
            h::FT
            ) where {FT<:AbstractFloat}
    return atmospheric_pressure(h) * FT(1e-6)
end








###############################################################################
#
# Calculate zenith angle
#
###############################################################################
"""
    zenith_angle(latd::FT, decd::FT, lhad::FT)
    zenith_angle(latd::FT, day::Int, hour::Int, minute::FT)

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
            day::Int,
            hour::Int,
            minute::FT=FT(0)
            ) where {FT<:AbstractFloat}
    _deg::FT = 360 / FT(YEAR_D) * (day + (hour+minute/60) / 24 + 10);
    decd::FT = -FT(23.44) * cosd(_deg);
    lhad::FT = (hour-12) * FT(15);

    return zenith_angle(latd, decd, lhad)
end




#end # module
