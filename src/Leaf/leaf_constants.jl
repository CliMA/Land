using Parameters


# Scaling functions for Photosynthesis temperature response and inhibition
ft(tl, ha) = exp(ha/(physcon.Rgas*(physcon.tfrz+25)) * (1-(physcon.tfrz+25)/tl));
fth(tl, hd, se, fc) = fc / (1 + exp((-hd+se*tl)/(physcon.Rgas*tl)));
fth25(hd, se) = 1.0 + exp( (-hd + se * (physcon.tfrz+25.)) / (physcon.Rgas * (physcon.tfrz+25.)) );

# Structure with all parameter temperature dependencies of a Leaf (largely based on Bonan's Photosynthesis model but exported as struct)
@with_kw mutable struct leaf_params
    # Rate constants (arbitrary units here, only relative rates are important):
    Kf = 0.05                           # Rate constant for fluorescence (might try to fit this eventually)
    Kd = 0.85                           # Rate constant for thermal dissipation
    Kp = 4.0                            # Rate constant for photochemistry (all reaction centers open)
    maxPSII = Kp/(Kp+Kf+Kd)             # Max PSII yield (Kn=0, all RC open)

    #use_acclim = false
    o₂ = 0.209e3                         # Standard O2 in mmol/mol
    T = 298.15                           # standard temperature (25C)

    # Curvature parameter:
    θ_j = 0.90

    # Conductances:
    # Ball Berry model
    go = 0.01
    g1 = 9

    # Add already here: Mesophyll conductance (high for now):
    gm = 1.0e99                             # Mesophyll conductance (μmol/m2/s/Pa)

    #---------------------------------------------------------------------
    # kc, ko, cp at 25C: Bernacchi et al (2001) Plant, Cell Environment 24:253-259
    # Derive sco from cp with o2=0.209 mol/mol and re-calculate Γ to allow
    # variation in O2
    #---------------------------------------------------------------------
    kc_25 = 404.9                         # Michaelis-Menten constant for CO2 at 25C μmol/mol
    ko_25 = 278.4                         # Michaelis-Menten constant for O2 at 25C  mmol/mol
    Γ_25_ = 42.75                         # CO2 compensation point at 25C μmol/mol
    sco = 0.5 * 0.209 / (Γ_25_ * 1.e-06)  # Γ_25 (μmol/mol) -> (mol/mol)
    Γ_25 = 0.5 * o₂ / sco * 1000.         # ! O2 is mmol/mol. Multiply by 1000 for μmol/mol

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
    kc  = kc_25
    ko  = ko_25
    Γ⋆  = Γ_25

    # Use some standard values first
    vcmax25 = 50                        # Leaf maximum carboxylation rate at 25C for canopy layer (μmol/m2/s)
    jmax25  = 1.6*vcmax25               # C3 - maximum electron transport rate at 25C for canopy layer (μmol/m2/s)
    rd25    = 0.7605                    # Leaf respiration rate at 25C for canopy layer (μmol CO2/m2/s)
    vcmax  = vcmax25
    jmax   = jmax25
    rdleaf = rd25
end

# Set Leaf rates with vcmax, jmax and rd at 25C as well as actual T here:
function setLeafT!(l::leaf_params,  T)
    l.T = T
    l.kc     = kc25          * ft(T, l.kcha)
    l.ko     = ko25          * ft(T, l.koha)
    l.Γ⋆     = cp25          * ft(T, l.cpha)
    l.vcmax   = l.vcmax25 * ft(T, l.vcmaxha) * fth(T, l.vcmaxhd, l.vcmaxse, l.vcmaxc)
    l.jmax    = l.jmax25  * ft(T, l.jmaxha)  * fth(T, l.jmaxhd, l.jmaxse, l.jmaxc)
    l.rdleaf  = l.rd25    * ft(T, l.rdha)    * fth(T, l.rdhd, l.rdse, l.rdc)
    # l.kd = max(0.8738,  0.0301*(T-273.15)+ 0.0773); # Can implement that later.
end
