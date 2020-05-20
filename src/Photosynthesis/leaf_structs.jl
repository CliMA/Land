
"""
    leaf_params{FT<:AbstractFloat}
    A structure, which contains all parameters for photosynthesis calculations of a single leaf

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct leaf_params{TT<:AbstractFloat}

    # broadband albedo and emissivity
    α::TT  = -999;                           # broadband shortwave albedo - absurd value to make sure we initialize correctly
    ε::TT  = -999;                           # longwave emissivity

    # thermal characteristics
    "DRY leaf mass area `[kg/m2]`"
    LMA::TT         = 100e-3;                # DRY leaf mass area, NOT MEAN VALUE WHICH INCLUDES RELATIVE WATER CONTENT - kg/m2 - density times thickness
    "leaf relative water content `[-]`"
    RWC::TT         = 0.8;                   # leaf relative water content
    "leaf specific heat capacity `[J/kg/K]`"
    Cleaf::TT       = 1000.0;                # leaf specific heat (J/kg/K)

    # pore traits
    #pore_density::TT   = 100.0  / (1e-6);     # pore density in pores/m^2

    "Chloroplast volume airspace `[m^3 per pore]`"
    Chloroplast_rel_volume::TT  = 2.0/100.0;        # hloroplast volume airspace (m^3 per pore) - http://kirschner.med.harvard.edu/files/bionumbers/A%20semi%20quantitative%20analysis%20of%20the%20photosynthetic%20system%20in%20an%20'average'%20C3%20plant%20leaf.pdf

    "Fraction of absorbed light used by PSII ETR"
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
    
    "dark-adapted fluorescence yield Fo (`Kp=max`)"
    Fo::TT = 0;     # dark-adapted fluorescence yield Fo,0
    "light-adapted fluorescence yield Fo' (`Kp=max`)"
    Fo′::TT = 0;     # light-adapted fluorescence yield in the dark Fo
    "dark adapted Fm yield (`Kp=0`)"
    Fm::TT = 0;     # light-adapted fluorescence yield Fm
    "light adapted Fm' yield (`Kp=0`)"
    Fm′::TT = 0;    # dark-adapted fluorescence yield Fm
    "steady-state (light-adapted) yield Ft (aka Fs)"
    ϕs::TT = 0;     # steady-state (light-adapted) yield Ft (aka Fs)

    "PSII yield"
    φ::TT = 0;

    "Photochemical quenching"
    qQ::TT = 0;     # 
    qE::TT = 0;     # non-photochemical quenching
    "Non-Photochemical quenching "
    NPQ::TT = 0;
    "Dynamic Steady state `[Bool]`"
    dynamic_state::Bool = false

    "Efficiency of photosynthesis light reactions`[mol CO₂/mol e-]`"
    effcon::TT = 1/5                       # [mol CO2/mol e-]  number of CO2 per electrons - typically 1/5 for C3 and 1/6 for C4, i.e. about 1/10 for both PSI and PSII!

    "Total efficiency, incl. photorespiration `[mol CO₂/mol e-]`"
    CO2_per_electron::TT = 1/5             # Similar to above but accounting for Photorespiration as well (which "steals" electrons)
    
    "Actual electron transport rate `[μmol m-2 s-1]`"
    Ja::TT = 0                               # Actual electron transport rate (needed to compute fluorescence later on)
    
    "O₂ mixing ratio `[mol/mol]`"
    o₂::TT = 0.209;                         # Standard O2 in mol/mol
    
    "Leaf temperature `[K)]`"
    T::TT = 298.15;                           # standard temperature (25C)
    
    "Previous Leaf temperature `[K]`"   # Will be checked, if no update, don't recalculate!
    T_old::TT = -1;
        
    "Saturation vapor pressure at leaf temperature `[Pa]`"
    esat::TT = 0.0

    "Derivative of saturation vapor pressure at leaf temperature wrt T `[Pa/K]`"
    desat::TT = 0.0
    #(esat, desat) = SatVap(T)
    "Curvature parameter:"
    θ_j::TT = 0.995;

    "total leaf conductance `[mmol H₂O /m2/s]`"
    gleaf::TT = 0.1;                           # total leaf conductance (μmol/m2/s)
    "Mesophyll conductance `[mmol H₂O /m2/s]`"
    gm::TT = Inf;                              # Mesophyll conductance (μmol/m2/s/)
    "Stomatal conductance `[mmol H₂O /m2/s]`"
    gs::TT = 0.1;                              # Stomatal conductance (μmol/m2/s); just using some prior
    "Steady state Stomatal conductance `[mmol H₂O /m2/s]`"
    gs_ss::TT = 0.1;                           # Steady state Stomatal conductance (μmol/m2/s);


    "Gross Photosynthetic Rate `[μmol/m2/s]`"
    Ag::TT = 0
    "Net Photosynthetic Rate `[μmol/m2/s]`"
    An::TT = 0
    "Rubisco Limited  Photosynthetic Rate `[μmol/m2/s]`"
    Ac::TT = 0
    "Light Limited  Photosynthetic Rate `[μmol/m2/s]`"
    Aj::TT = 0
    "Product Limited  Photosynthetic Rate `[μmol/m2/s]`"
    Ap::TT = 0

    "Sensible Heat Flux `[W/m2]`"
    H::TT = 0
    "Latent Heat Flux `[W/m2]`"
    LE::TT = 0
    "Net Radiation Balance `[W/m2]`"
    Rn::TT = 0

    "Sap flow `[?]`"
    Sap::TT = 0;
    "Absorbed Photosynthetically Active Radiation `[μmol/m2/s]`"
    APAR::TT = 100
             
    # Placeholders for MM constants (need to be set with current Temperature):
    "Michaelis-Menten constant for CO₂ `[Pa]`"
    Kc::TT  = 0;
    "Michaelis-Menten constant for O₂ `[Pa]`"
    Ko::TT  = 0;
    "Michaelis-Menten constant for PEP carboxylase `[Pa]`"
    Kpep::TT  = 0;
    "CO₂ compensation point `[Pa]`"
    Γstar::TT  = 0;

    # Use some standard values first
    "Leaf maximum carboxylation rate at 25C `[μmol/m2/s]`"
    Vcmax25::TT = 80.0;                        
    "Leaf maximum PEP carboxylation rate at 25C `[μmol/m2/s]`"
    Vpmax25::TT = 120.0;                        
    "Maximum electron transport rate at 25C `[μmol/m2/s]`"
    Jmax25::TT  = 1.97*Vcmax25;                
    "Leaf respiration rate at 25C `[μmol/m2/s]`"
    Rd25::TT    = 0.01*Vcmax25;                      
    
    "Actual Leaf maximum carboxylation rate `[μmol/m2/s]`"
    Vcmax::TT  = Vcmax25;
    "Actual Leaf maximum PEP carboxylase rate `[μmol/m2/s]`"
    Vpmax::TT  = Vpmax25;
    "Actual Maximum electron transport rate `[μmol/m2/s]`"
    Jmax::TT   = Jmax25;
    "Actual Leaf respiration rate `[μmol/m2/s]`"
    Rd::TT = Rd25;

    "CO₂ concentration at leaf surface `[μmol/mol]`"
    Cs::TT = 400.0;  
    "CO₂ concentration in chloroplast `[μmol/mol]"
    Cc::TT = 400.0;                                                   
    "relative humidity at the surface of the leaf `[-]`"
    RH::TT = 100.0;                               
    "Vapor pressure difference across leaf interface `[Pa]`"
    VPD::TT = 0.1;                                # VPD gradient ACROSS THE LEAF interface - NOT the weather station one 8-)

    # to be computed
    "Electron transport rate `[μmol/m2/s]`, light limited only"
    Je::TT = 0.0 ;                              # electron transport rate

    "Potential Electron Transport Rate `[μmol/m2/s]`"
    Je_pot::TT = 0.0
    
    # plant hydraulics
    "leaf water potential `[Pa]`"
    psi_l::TT   = -1.5e6;                         # leaf water potential (Pa)
    "leaf water potential at 50% drop in conductivity `[Pa]`"
    psi_l50::TT = -1.75e6;                        # leaf water potential at 50% drop in conductivity (Pa)
    "maximum leaf-xylem conductivity `[m/s]`"
    kmax::TT    = 4e-8                            # maximum leaf-xylem conductivity (m/s)
    "actual xylem conductivity `[m/s]`"
    kx::TT      = kmax                            # actual xylem conductivity (m/s)
    "leaf level hydraulic conductivity `[m/s]`"
    kleaf::TT   = kmax                            # leaf level hydraulic conductivity (m/s)
    "slope of Weibull curve `[m/s/Pa]`"
    ck::TT      = 2.95;                           # slope of Weibull curve
    #ε_modulus::TT = 20e6;                        # elastic modulus - for later

    "tree height `[m]`"
    height      = 20.;                            # 
    "tree capacitance `[kg/m3/Pa]`"
    Ctree::TT   = (79 + 8*height)*1e-6;           # tree capacitance (kg m−3 Pa−1)
                                                  # -Scholz et al. 2011 Book, Hydraulic Capacitance: Biophysics and Functional Significance,
                                                  # can also be related to P50 see same book

    
    use_colim::Bool = true                    # Use co-limitation

    
    z0m         = -999.;                          # tree roughness (m)
    z0h         = -999.;                          # tree roughness (m) - TODO should be changed later
    d           = -999.;                          # tree displacement height (m)

    dleaf       = 2e-3;                           # leaf thickness (m)
    LAI         = 1.0;                            # leaf area thickness (m^2/m^2)
    Cd          = 0.01;                           # m/sqrt(s) turbulent transfer coefficient
end

"""
    meteo
    A Placeholder for all Meteo data for the moment

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct meteo{TT<:Number}

    S_down::TT = 0.;
    L_down::TT = 0.;
    "Air temperature `[K]`"
    T_air::TT  = 280.;      # T in K
    "Air vapor pressure `[Pa]`"
    e_air::TT  = 1500.;
    "Air surface pressure `[Pa]`"
    P_air::TT  =  1e5 ; 
    "Ambient CO₂ concentration `[µmol/mol]`"
    Ca::TT     =  400.;
    "PAR `[µmol/m2/s]`"
    PAR::TT    = 0.;
    "Wind Speed `[m/s]`"
    U::TT      = 1e-6;
    "Measurement height `[m]`"
    zscreen::TT= 10.0; # measurement height - default
    "Measurement height `[m]`"
    L::TT      = 1e6;  # atmospheric Obukhov length
    # parameter to define stability function for stable case
    stab_type_stable::TT = 1; # 2 Webb correction tends to be unstable at night - suggest not using
    ustar::TT = 1e-6
    "Conversion factor from `[m/s]` to `[mol/m2/s]`"
    g_m_s_to_mol_m2_s::TT = -Inf
    "Conversion factor from mixing ratio `[μmol/mol]` to partial pressure `[Pa]`"
    ppm_to_Pa::TT = 0.1
    "Air resistance `[Pa]`"
    ra::TT    = 1e6
    "absorbed PAR `[µmol/m2/s]`"
    APAR::TT = 500.0
    "CO₂ concentration at leaf surface `[µmol/mol]`"
    Cs::TT = 0.0
end


# Need to ba abstracted as well!
function setkx!(l::leaf_params, psis, psi_l) # set hydraulic conductivitytimes Delta Psi
    l.kx = l.kmax * IntWeibull(psis,psi_l,l.psi_l50,l.ck)/max(psis-psi_l,1e-6); # kmax . int_psis^psil k(x)dx = kmax . IntWeibull(psil);
    #println("k_xylem = ",l.kx," psi_s=",psis," psi_l=",psi_l)
end

function setLeafkl!(l::leaf_params, psi_l) # set hydraulic conductivity
    l.kleaf = l.kmax * Weibull(psi_l,l.psi_l50,l.ck); # kmax . int_psis^psil k(x)dx = kmax . IntWeibull(psil);
end

