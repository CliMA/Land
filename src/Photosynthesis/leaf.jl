
# 
"""
    leaf_params
    Structure with all parameter temperature dependencies of a Leaf (largely based on Bonan's Photosynthesis model but exported as struct)

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct leaf_params{TT<:Number}

    # broadband albedo and emissivity
    α::TT  = -999;                           # broadband shortwave albedo - absurd value to make sure we initialize correctly
    ε::TT  = -999;                           # longwave emissivity

    # thermal characteristics
    "DRY leaf mass area"
    LMA::TT         = 100e-3;                # DRY leaf mass area, NOT MEAN VALUE WHICH INCLUDES RELATIVE WATER CONTENT - kg/m2 - density times thickness
    "leaf relative water content"
    RWC::TT         = 0.8;                   # leaf relative water content
    "leaf specific heat (J/kg/K)"
    Cleaf::TT       = 1000.0;                # leaf specific heat (J/kg/K)

    # pore traits
    #pore_density::TT   = 100.0  / (1e-6);     # pore density in pores/m^2
    Chloroplast_rel_volume::TT  = 2.0/100.0;        # hloroplast volume airspace (m^3 per pore) - http://kirschner.med.harvard.edu/files/bionumbers/A%20semi%20quantitative%20analysis%20of%20the%20photosynthetic%20system%20in%20an%20'average'%20C3%20plant%20leaf.pdf

    "Fraction of  absorbed light used by PSII ETR"
    PSII_frac::TT = 0.5
    # Rate constants (arbitrary units here, only relative rates are important):
    "Rate constant for fluorescence (const)"
    Kf::TT = 0.05;                           # Rate constant for fluorescence (might try to fit this eventually)
    "Rate constant for thermal dissipation"
    Kd::TT = 0.85;                           # Rate constant for thermal dissipation, can include T dependence later!
    "Rate constant for photochemistry (all reaction centers open)"
    Kp::TT = 4.0;                            # Rate constant for photochemistry (all reaction centers open)
    "NPQ rate constant"
    Kn::TT = 0;                              # NPQ rate constant (initially zero)
    Kn_ss::TT = 0;                           # Kn in steady state
    "max PSII yield"
    maxPSII::TT = Kp/(Kp+Kf+Kd);             # Max PSII yield (Kn=0, all RC open)
    "Flexas et al derived Fluorescence model params"
    Knparams = [5.01, 1.93, 10.0]
    "dark-adapted fluorescence yield Fo"
    Fo::TT = 0;     # dark-adapted fluorescence yield Fo,0
    "light-adapted fluorescence yield in the dark Fo'"
    Fo′::TT = 0;     # light-adapted fluorescence yield in the dark Fo
    "dark adapted Fm yield"
    Fm::TT = 0;     # light-adapted fluorescence yield Fm
    "light adapted Fm' yield"
    Fm′::TT = 0;    # dark-adapted fluorescence yield Fm
    "steady-state (light-adapted) yield Ft (aka Fs)"
    ϕs::TT = 0;     # steady-state (light-adapted) yield Ft (aka Fs)
    "fluorescence efficiency (to be used for SIF RT)"
    eta::TT = 0;
    qQ::TT = 0;     # photochemical quenching
    qE::TT = 0;     # non-photochemical quenching
    NPQ::TT = 0;
    "Dynamic Steady state boolean"
    dynamic_state::Bool = false


    effcon::TT = 1/5                       # [mol CO2/mol e-]  number of CO2 per electrons - typically 1/5 for C3 and 1/6 for C4, i.e. about 1/10 for both PSI and PSII!
    CO2_per_electron::TT = 1/5             # Similar to above but accounting for Photorespiration as well (which "steals" electrons)
    "Actual electron transport rate (μmol m-2 s-1)"
    Ja::TT = 0                               # Actual electron transport rate (needed to compute fluorescence later on)
    "C3 or C4 plant boolean"
    C3::Bool = true                           # C3 (or C4) plant
    use_colim::Bool = true                    # Use co-limitation

    gstyp::Int = 1                            # stomatal conductance type (1=Ball Berry, 0 = Medlyn)
    "O2 (mol/mol)"
    o₂::TT = 0.209;                         # Standard O2 in mol/mol
    "Leaf temperature (K)"
    T::TT = 298.15;                           # standard temperature (25C)
    esat::TT = 0.0
    desat::TT = 0.0
    #(esat, desat) = SatVap(T)
    # Curvature parameter:
    θ_j::TT = 0.90;

    # Conductances:
    # Ball Berry or Medlyn model
    g0::TT = 0.1;                              # Ball-Berry minimum leaf conductance, unstressed (mol/m2/s) - Sellers  et  al.  (1996)  used  0.1 b=  for  C3  plants  and  0.4 b=for  C4  plants
    g1_BB::TT = 9.0;                           # Ball-Berry slope of conductance-photosynthesis relationship, unstressed - m=  for  C3  plants  and  4m=  for  C4  plants  (Collatz  et  al.  1991,  1992)
    g1_Medlyn::TT = 126.49;                    # Medlyn slope of conductance-photosynthesis relationship, unstressed - Pa^(1/2) Medlyn et al. 2017
    vpd_min = 100.0                            # Min VPD for Medlyn model (blows up otherwise?)
    # Add already here: Mesophyll conductance (Inf for now):

    "total leaf conductance (μmol/m2/s)"
    gleaf::TT = 0.1;                           # total leaf conductance (μmol/m2/s)
    "Mesophyll conductance (μmol/m2/s/)"
    gm::TT = Inf;                              # Mesophyll conductance (μmol/m2/s/)
    "Stomatal conductance (μmol/m2/s)"
    gs::TT = 0.1;                              # Stomatal conductance (μmol/m2/s); just using some prior
    "Steady state Stomatal conductance (μmol/m2/s)"
    gs_ss::TT = 0.1;                           # Steady state Stomatal conductance (μmol/m2/s);

             
    # Placeholders for MM constants (need to be set with current Temperature):
    "Michaelis-Menten constant for CO₂ (Pa)"
    Kc::TT  = 0;
    "Michaelis-Menten constant for O₂ (Pa)"
    Ko::TT  = 0;
    "Michaelis-Menten constant for PEP carboxylase (Pa)"
    Kpep::TT  = 0;
    "CO₂ compensation point (Pa)"
    Γstar::TT  = 0;

    # Use some standard values first
    "Leaf maximum carboxylation rate at 25C (μmol/m2/s)"
    Vcmax25::TT = 80.0;                        # Leaf maximum carboxylation rate at 25C for canopy layer (μmol/m2/s)
    "Leaf maximum PEP carboxylation rate at 25C (μmol/m2/s)"
    Vpmax25::TT = 120.0;                        # Leaf maximum carboxylation rate at 25C for canopy layer (μmol/m2/s)
    "Maximum electron transport rate at 25C (μmol/m2/s)"
    Jmax25::TT  = 1.97*Vcmax25;                # C3 - maximum electron transport rate at 25C for canopy layer (μmol/m2/s)
    "Leaf respiration rate at 25C (μmol CO2/m2/s)"
    Rd25::TT    = 0.01*Vcmax25;                      # Leaf respiration rate at 25C for canopy layer (μmol CO2/m2/s)
    
    "Actual Leaf maximum carboxylation rate (μmol/m2/s)"
    Vcmax::TT  = Vcmax25;
    "Actual Leaf maximum PEP carboxylase rate (μmol/m2/s)"
    Vpmax::TT  = Vpmax25;
    "Actual Maximum electron transport rate (μmol/m2/s)"
    Jmax::TT   = Jmax25;
    "Actual Leaf respiration rate (μmol CO2/m2/s)"
    Rdleaf::TT = Rd25;

    "CO2 concentration in chloroplast (ppm)"
    Cc::TT = 400.0;                               
    "CO2 concentration in chloroplast (Pa)"
    Cc_Pa::TT = 40;                               
    "relative humidity at the surface of the leaf"
    RH::TT = 100.0;                               
    "Vapor pressure difference across leaf interface (kPa)"
    VPD::TT = 0.1;                                # VPD gradient ACROSS THE LEAF interface - NOT the weather station one 8-)

    # to be computed
    "Electron transport rate (μmol m-2 s-1)"
    Je::TT = 0.0 ;                              # electron transport rate
    # tree/leaf traits
    height      = 20.;                            # tree height (m)
    z0m         = -999.;                          # tree roughness (m)
    z0h         = -999.;                          # tree roughness (m) - TODO should be changed later
    d           = -999.;                          # tree displacement height (m)

    dleaf       = 2e-3;                           # leaf thickness (m)
    LAI         = 1.0;                            # leaf area thickness (m^2/m^2)
    Cd          = 0.01;                           # m/sqrt(s) turbulent transfer coefficient

    # plant hydraulics
    psi_l::TT   = -1.5e6;                         # leaf water potential (Pa)
    psi_l50::TT = -1.75e6;                        # leaf water potential at 50% drop in conductivity (Pa)
    kmax::TT    = 4e-8                            # maximum leaf-xylem conductivity (m/s)
    kx::TT      = kmax                            # actual xylem conductivity (m/s)
    kleaf::TT   = kmax                            # leaf level hydraulic conductivity (m/s)
    ck::TT      = 2.95;                           # slope of Weibull curve
    #ε_modulus::TT = 20e6;                        # elastic modulus - for later
    Ctree::TT   = (79 + 8*height)*1e-6;           # tree capacitance (kg m−3 Pa−1)
                                                  # -Scholz et al. 2011 Book, Hydraulic Capacitance: Biophysics and Functional Significance,
                                                  # can also be related to P50 see same book
end

# Set Leaf rates with Vcmax, Jmax and rd at 25C as well as actual T here:
# For some reason, this is slow and allocates a lot, can be improved!!
"Set Leaf rates with Vcmax, Jmax and rd at 25C as well as actual T here"
function setLeafT!(l::leaf_params)
    l.Kc      = l.Kc_25          * ft(l.T, l.Kcha);
    l.Ko      = l.Ko_25          * ft(l.T, l.Koha);
    l.Γstar   = l.Γ_25           * ft(l.T, l.cpha);
    l.Vcmax   = l.Vcmax25 * ft(l.T, l.Vcmaxha) * fth(l.T, l.Vcmaxhd, l.Vcmaxse, l.Vcmaxc);
    l.Jmax    = l.Jmax25  * ft(l.T, l.Jmaxha)  * fth(l.T, l.Jmaxhd, l.Jmaxse, l.Jmaxc);
    l.Rdleaf  = l.Rd25    * ft(l.T, l.rdha)    * fth(l.T, l.rdhd, l.rdse, l.rdc);
    (l.esat, l.desat) = SatVap(l.T);
    # l.kd = max(0.8738,  0.0301*(l.T-273.15)+ 0.0773); # Can implement that later.
end

function setkx!(l::leaf_params, psis, psi_l) # set hydraulic conductivitytimes Delta Psi
    l.kx = l.kmax * IntWeibull(psis,psi_l,l.psi_l50,l.ck)/max(psis-psi_l,1e-6); # kmax . int_psis^psil k(x)dx = kmax . IntWeibull(psil);
    #println("k_xylem = ",l.kx," psi_s=",psis," psi_l=",psi_l)
end

function setLeafkl!(l::leaf_params, psi_l) # set hydraulic conductivity
    l.kleaf = l.kmax * Weibull(psi_l,l.psi_l50,l.ck); # kmax . int_psis^psil k(x)dx = kmax . IntWeibull(psil);
end

