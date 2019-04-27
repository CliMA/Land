using Parameters

@with_kw struct phys
    grav = 9.80665;               # Gravitational acceleration (m/s2)
    tfrz = 273.15;                # Freezing point of water (K)
    σ = 5.67e-08;                 # Stefan-Boltzmann constant (W/m2/K4)
    mmdry = 28.97 / 1000;         # Molecular mass of dry air (kg/mol)
    mmh2o = 18.02 / 1000;         # Molecular mass of water (kg/mol)
    Cpd = 1005;                   # Specific heat of dry air at constant pressure (J/kg/K)
    Cpw = 1846;                   # Specific heat of water vapor at constant pressure (J/kg/K)
    Rgas = 8.31446;               # Universal gas constant (J/K/mol)
    visc_0 = 13.3e-06;            # Kinematic viscosity at 0C and 1013.25 hPa (m2/s)
    Dh_0 = 18.9e-06;              # Molecular diffusivity (heat) at 0C and 1013.25 hPa (m2/s)
    Dv_0 = 21.8e-06;              # Molecular diffusivity (H2O) at 0C and 1013.25 hPa (m2/s)
    Dc_0 = 13.8e-06;              # Molecular diffusivity (CO2) at 0C and 1013.25 hPa (m2/s)
end

physcon = phys()

# Scaling functions for Photosynthesis temperature response and inhibition
ft(tl, ha) = exp(ha/(physcon.Rgas*(physcon.tfrz+25)) * (1-(physcon.tfrz+25)/tl));
fth(tl, hd, se, fc) = fc / (1 + exp((-hd+se*tl)/(physcon.Rgas*tl)));
fth25(hd, se) = 1.0 + exp( (-hd + se * (physcon.tfrz+25.)) / (physcon.Rgas * (physcon.tfrz+25.)) );

# Structure with all parameter temperature dependencies of a Leaf (largely based on Bonan's Photosynthesis model but exported as struct)
@with_kw mutable struct leaf_params
    #use_acclim = false
    o2 = 0.209e3                         # Standard O2 in mmol/mol
    T = 298.15                           # standard temperature (25C)
    #---------------------------------------------------------------------
    # kc, ko, cp at 25C: Bernacchi et al (2001) Plant, Cell Environment 24:253-259
    # Derive sco from cp with o2=0.209 mol/mol and re-calculate Γ to allow
    # variation in O2
    #---------------------------------------------------------------------
    kc_25 = 404.9                         # Michaelis-Menten constant for CO2 at 25C μmol/mol
    ko_25 = 278.4                         # Michaelis-Menten constant for O2 at 25C  mmol/mol
    Γ_25_ = 42.75                         # CO2 compensation point at 25C μmol/mol
    sco = 0.5 * 0.209 / (Γ_25_ * 1.e-06)  # Γ_25 (μmol/mol) -> (mol/mol)
    Γ_25 = 0.5 * o2 / sco * 1000.         # ! O2 is mmol/mol. Multiply by 1000 for μmol/mol

    #---------------------------------------------------------------------
    # Activation energy:
    # Bernacchi et al (2001) Plant, Cell Environment 24:253-259
    # Bernacchi et al (2003) Plant, Cell Environment 26:1419-1430
    # Acclimation from: Kattge and Knorr (2007) Plant, Cell Environment 30:1176-1190
    #---------------------------------------------------------------------
    kcha    = 79430.                      # Activation energy for kc (J/mol)
    koha    = 36380.                      # Activation energy for ko (J/mol)
    cpha    = 37830.                      # Activation energy for cp (J/mol)
    vcmaxha = 65330.                      # Activation energy for Vcmax (J/mol)
    jmaxha  = 43540.                      # Activation energy for Jmax (J/mol)
    rdha    = 46390.                      # Activation energy for Rd (J/mol)

    #---------------------------------------------------------------------
    # High temperature deactivation:
    # Leuning (2002) Plant, Cell Environment 25:1205-1210
    # The factor "c" scales the deactivation to a value of 1.0 at 25C
    # Acclimation from: Kattge and Knorr (2007) Plant, Cell Environment 30:1176-1190
    #---------------------------------------------------------------------
    vcmaxhd = 150000.                    # Deactivation energy for Vcmax (J/mol)
    jmaxhd  = 150000.                    # Deactivation energy for Jmax (J/mol)
    rdhd    = 150000.                    # Deactivation energy for Rd (J/mol)

    vcmaxse = 490.                       # Entropy term for Vcmax (J/mol/K)
    jmaxse  = 490.                       # Entropy term for Jmax (J/mol/K)
    rdse    = 490.                       # Entropy term for Rd (J/mol/K)

    vcmaxc = fth25(vcmaxhd, vcmaxse)    # Scaling factor for high temperature inhibition (25 C = 1.0)
    jmaxc  = fth25(jmaxhd, jmaxse)      # Scaling factor for high temperature inhibition (25 C = 1.0)
    rdc    = fth25(rdhd, rdse)          # Scaling factor for high temperature inhibition (25 C = 1.0)

    # Initialize at standard T:
    kc = kc_25
    ko = ko_25
    Γ  = Γ_25

    vcmax25 = 50                        # Leaf maximum carboxylation rate at 25C for canopy layer (μmol/m2/s)
    jmax25  = 1.6*vcmax25               # C3 - maximum electron transport rate at 25C for canopy layer (μmol/m2/s)
    rd25    = 0.7605                    # Leaf respiration rate at 25C for canopy layer (μmol CO2/m2/s)
    vcmax  = vcmax25
    jmax   = jmax25
    rdleaf = rd25
end

# Set Leaf rates with vcmax, jmax and rd at 25C as well as actual T here:
function setC3LeafRates!(l::leaf_params, vcmax25, jmax25, rd25, T)
    l.vcmax25 = vcmax25
    l.jmax25  = jmax25
    l.rd25    = rd25
    l.vcmax   = vcmax25 * ft(T, l.vcmaxha) * fth(T, l.vcmaxhd, l.vcmaxse, l.vcmaxc)
    l.jmax    = jmax25  * ft(T, l.jmaxha)  * fth(T, l.jmaxhd, l.jmaxse, l.jmaxc)
    l.rdleaf  = rd25    * ft(T, l.rdha)    * fth(T, l.rdhd, l.rdse, l.rdc)
end
