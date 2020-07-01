module Earth

using ..LandParameters

const GAS_R   = LandParameters.GAS_R
const GRAVITY = LandParameters.GRAVITY
const K_25    = LandParameters.K_25
const P_ATM   = LandParameters.P_ATM
const RK_25   = LandParameters.RK_25
const YEAR_D  = LandParameters.YEAR_D

# export public functions
export atmospheric_pressure,
       ppm_to_Pa,
       zenith_angle




###############################################################################
#
# Atmospheric pressure and partial pressure
#
###############################################################################
"""
    atmospheric_pressure(h::FT)

Calculate the atmospheric pressure, given
- `h` elevation in `[m]`
"""
function atmospheric_pressure(
            h::FT
            ) where {FT<:AbstractFloat}
    m::FT = 0.02896968;
    p::FT = FT(P_ATM) * exp( -FT(GRAVITY) * h * m / FT(RK_25) );

    return p
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
    zenith_angle(latd::FT, day::Number, hour::Number, minute::FT)

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
            day::Number,
            hour::Number,
            minute::FT=FT(0)
            ) where {FT<:AbstractFloat}
    _deg::FT = 360/FT(YEAR_D) * (day + (hour+minute/60) / 24 + 10);
    decd::FT = -FT(23.44) * cosd(_deg);
    lhad::FT = (hour-12) * FT(15);

    return zenith_angle(latd, decd, lhad)
end




end # module
