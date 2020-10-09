# TODO Move the constants to CLIMAParameters or LandParameters
Base.@kwdef struct phys{T<:Number}
    visc_0::T = 13.3e-06;            # Kinematic viscosity at 0C and 1013.25 hPa (m2/s)
    Dh_0::T = 18.9e-06;              # Molecular diffusivity (heat) at 0C and 1013.25 hPa (m2/s)
    Dv_0::T = 21.8e-06;              # Molecular diffusivity (H2O) at 0C and 1013.25 hPa (m2/s)
    Dc_0::T = 13.8e-06;              # Molecular diffusivity (CO2) at 0C and 1013.25 hPa (m2/s)
    Wtoμmole_s::T = 4.57;            # converts Watts/m2 to μmole/s
    Cdryleaf::T = 1396;              # Specific heat of dry leaf at constant pressure (J/kg/K)
end
