module PhysCon

using Parameters
#-----------------------------------------------------------------------
# DESCRIPTION:
# Calculate saturation vapor pressure and latent heat of vaporization
#
export physcon
Base.@kwdef struct phys{T<:Number}
    grav::T = 9.80665;               # Gravitational acceleration (m/s2)
    tfrz::T = 273.15;                # Freezing point of water (K)
    σ::T = 5.67e-08;                 # Stefan-Boltzmann constant (W/m2/K4)
    mmdry::T = 28.97 / 1000;         # Molecular mass of dry air (kg/mol)
    mmh2o::T = 18.02 / 1000;         # Molecular mass of water (kg/mol)
    Cpd::T = 1005.;                  # Specific heat of dry air at constant pressure (J/kg/K)
    Cpw::T = 1846.;                  # Specific heat of water vapor at constant pressure (J/kg/K)
    Cpl::T = 4188.;                  # Specific heat of liquid water  at constant pressure at 15°C (J/kg/K)
    Rgas::T = 8.31446;               # Universal gas constant (J/K/mol)
    visc_0::T = 13.3e-06;            # Kinematic viscosity at 0C and 1013.25 hPa (m2/s)
    Dh_0::T = 18.9e-06;              # Molecular diffusivity (heat) at 0C and 1013.25 hPa (m2/s)
    Dv_0::T = 21.8e-06;              # Molecular diffusivity (H2O) at 0C and 1013.25 hPa (m2/s)
    Dc_0::T = 13.8e-06;              # Molecular diffusivity (CO2) at 0C and 1013.25 hPa (m2/s)
    Rd::T   = Rgas/mmdry;            # ideal gas constant for dry air (J/kg/K)
    ε::T    = mmh2o/mmdry            # ratio of Molecular mass of water over Molecular mass of dry air
    ρw::T   = 1000.;                 # density of liquid water
    Wtoμmole_s::T = 4.57;            # converts Watts/m2 to μmole/s
    K::T    = 0.40;                  # von Karman constant
    Cdryleaf::T = 1396;              # Specific heat of dry leaf at constant pressure (J/kg/K)
end

physcon = phys{Float32}()

end
