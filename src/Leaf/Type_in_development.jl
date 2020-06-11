###############################################################################
#
# These types should go to type.jl when documentations are done
#
###############################################################################








###############################################################################
#
# Leaf photosynthesis-related parameter set
# This struct passed the FT test
# More documentation required
#
###############################################################################
"""
    LeafParams{FT<:AbstractFloat}

A structure that contains all parameters for photosynthesis calculations of a single leaf.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct LeafParams{FT<:AbstractFloat}
    # TODO Group these variables

    # broadband albedo and emissivity
    "broadband shortwave albedo"
    α::FT = -999
    "longwave emissivity"
    ε::FT = -999

    # thermal-related characteristics
    "Leaf dry mass per area `[kg m⁻²]`, not WET mass"
    LMA::FT = 0.1
    "leaf relative water content"
    RWC::FT = 0.8
    "leaf specific heat capacity `[J kg⁻¹ K⁻¹]`"
    Cleaf::FT = 1000.0

    # TODO remove? pore density in pores/m^2
    #pore_density::FT = 1e8

    # chloroplast volume airspace (m^3 per pore)
    # http://kirschner.med.harvard.edu/files/bionumbers/A%20semi%20quantitative%20analysis%20of%20the%20photosynthetic%20system%20in%20an%20'average'%20C3%20plant%20leaf.pdf
    "Chloroplast volume airspace `[m³ per pore]`"
    Chloroplast_rel_volume::FT = 2.0/100.0

    # Photosynthesis
    "Fraction of absorbed light used by PSII ETR"
    PSII_frac::FT = 0.5
    # Rate constants (arbitrary units here, only relative rates are important):
    # TODO try to fit this eventually?
    "Rate constant for fluorescence (const)"
    Kf::FT = 0.05
    # TODO include T dependence later?
    "Rate constant for thermal dissipation"
    Kd::FT = 0.85
    "Rate constant for photochemistry (all reaction centers open)"
    Kp::FT = 4.0
    "NPQ rate constant (initially zero)"
    Kn::FT = 0
    "Kn in steady state"
    Kn_ss::FT = 0
    "max PSII yield (Kn=0, all RC open)"
    maxPSII::FT = Kp/(Kp+Kf+Kd)

    # Fluorescence-related
    "dark-adapted fluorescence yield (`Kp=max`)"
    Fo::FT = 0
    "light-adapted fluorescence yield in the dark (`Kp=max`)"
    Fo′::FT = 0
    "dark adapted yield (`Kp=0`)"
    Fm::FT = 0
    "light adapted yield (`Kp=0`)"
    Fm′::FT = 0
    "steady-state (light-adapted) yield (aka Fs)"
    ϕs::FT = 0

    # TODO Is this the quntum yield of electron?
    "PSII yield"
    φ::FT = 0;

    "Photochemical quenching"
    qQ::FT = 0
    "energy quenching"
    qE::FT = 0;
    "Non-Photochemical quenching "
    NPQ::FT = 0;
    "Dynamic Steady state `[Bool]`"
    dynamic_state::Bool = false

    # number of CO₂ per electrons - typically 1/5 for C3 and 1/6 for C4, i.e. about 1/10 for both PSI and PSII!
    "Efficiency of photosynthesis light reactions`[mol CO₂/mol e-]`"
    effcon::FT = 1/5

    # Similar to `effcon` but accounting for Photorespiration as well (which "steals" electrons)
    "Total efficiency, incl. photorespiration `[mol CO₂/mol e-]`"
    CO₂_per_electron::FT = 1/5

    # Actual electron transport rate (needed to compute fluorescence later on)
    "Actual electron transport rate `[μmol m-2 s-1]`"
    Ja::FT = 0

    "O₂ mixing ratio `[mol/mol]`"
    O₂::FT = 0.209
    "O₂ partial pressure `[Pa]`"
    p_O₂::FT = 21176.926

    "Leaf temperature `[K)]`"
    T::FT = 298.15

    # Will be checked, if no update, don't recalculate!
    "Previous Leaf temperature `[K]`"
    T_old::FT = -1;

    "Saturation vapor pressure at leaf temperature `[Pa]`"
    esat::FT = 0.0

    "Derivative of saturation vapor pressure at leaf temperature wrt T `[Pa/K]`"
    desat::FT = 0.0

    "Curvature parameter:"
    θ_j::FT = 0.995;

    "total leaf diffusive conductance `[mol H₂O m⁻² s⁻¹]`"
    gleaf::FT = 0.1
    "Mesophyll conductance `[mol H₂O m⁻² s⁻¹]`"
    gm::FT = Inf
    "Stomatal conductance `[mol H₂O m⁻² s⁻¹]`"
    gs::FT = 0.1
    "Steady state Stomatal conductance `[mol H₂O m⁻² s⁻¹]`"
    gs_ss::FT = 0.1

    "Gross Photosynthetic Rate `[μmol m⁻² s⁻¹]`"
    Ag::FT = 0
    "Net Photosynthetic Rate `[μmol m⁻² s⁻¹]`"
    An::FT = 0
    "Rubisco Limited Photosynthetic Rate `[μmol m⁻² s⁻¹]`"
    Ac::FT = 0
    "Light Limited Photosynthetic Rate `[μmol m⁻² s⁻¹]`"
    Aj::FT = 0
    "Product Limited Photosynthetic Rate `[μmol m⁻² s⁻¹]`"
    Ap::FT = 0

    "Sensible Heat Flux `[W m⁻²]`"
    H::FT = 0
    "Latent Heat Flux `[W m⁻²]`"
    LE::FT = 0
    "Net Radiation Balance `[W m⁻²]`"
    Rn::FT = 0

    # TODO Update this with Plant Hydraulics
    "Sap flow `[?]`"
    Sap::FT = 0;
    # TODO Remove this if not necessary? Or update it with time?
    "Absorbed Photosynthetically Active Radiation `[μmol m⁻² s⁻¹]`"
    APAR::FT = 100

    # Placeholders for MM constants (need to be set with current Temperature):
    "Michaelis-Menten constant for CO₂ `[Pa]`"
    Kc::FT = 0;
    "Michaelis-Menten constant for O₂ `[Pa]`"
    Ko::FT = 0;
    "Michaelis-Menten constant for PEP carboxylase `[Pa]`"
    Kpep::FT = 0;
    "CO₂ compensation point `[Pa]`"
    Γstar::FT = 0;

    # Use some standard values first
    "Leaf maximum carboxylation rate at 25C `[μmol m⁻² s⁻¹]`"
    Vcmax25::FT = 80.0
    "Leaf maximum PEP carboxylation rate at 25C `[μmol m⁻² s⁻¹]`"
    Vpmax25::FT = 120.0
    "Maximum electron transport rate at 25C `[μmol m⁻² s⁻¹]`"
    Jmax25::FT = 1.97*Vcmax25
    "Leaf respiration rate at 25C `[μmol m⁻² s⁻¹]`"
    Rd25::FT = 0.01*Vcmax25

    "Actual Leaf maximum carboxylation rate `[μmol m⁻² s⁻¹]`"
    Vcmax::FT = Vcmax25;
    "Actual Leaf maximum PEP carboxylase rate `[μmol m⁻² s⁻¹]`"
    Vpmax::FT = Vpmax25;
    "Actual Maximum electron transport rate `[μmol m⁻² s⁻¹]`"
    Jmax::FT = Jmax25;
    "Actual Leaf respiration rate `[μmol m⁻² s⁻¹]`"
    Rd::FT = Rd25;

    "CO₂ partial pressure at leaf surface `[Pa]`"
    p_s::FT = 40.0
    "leaf internal CO₂ partial pressure `[Pa]`"
    p_i::FT = 30.0

    "CO₂ concentration at leaf surface `[μmol/mol]`"
    Cs::FT = 400.0
    "CO₂ concentration in chloroplast `[μmol/mol]"
    Cc::FT = 400.0
    "relative humidity at the surface of the leaf `[-]`"
    RH::FT = 100.0

    # VPD gradient ACROSS THE LEAF interface - NOT the weather station one 8-)
    "Vapor pressure difference across leaf interface `[Pa]`"
    VPD::FT = 0.1

    # to be computed
    "Electron transport rate `[μmol m⁻² s⁻¹]`, light limited only"
    Je::FT = 0.0
    "Potential Electron Transport Rate `[μmol m⁻² s⁻¹]`"
    Je_pot::FT = 0.0

    # TODO May need to remove to avoid duplication
    # plant hydraulics
    "leaf water potential `[Pa]`"
    psi_l::FT = -1.5e6
    "leaf water potential at 50% drop in conductivity `[Pa]`"
    psi_l50::FT = -1.75e6
    "maximum leaf-xylem conductivity `[m s⁻¹]`"
    kmax::FT = 4e-8
    "actual xylem conductivity `[m s⁻¹]`"
    kx::FT = kmax
    "leaf level hydraulic conductivity `[m s⁻¹]`"
    kleaf::FT = kmax
    "slope of Weibull curve `[m s⁻¹/Pa]`"
    ck::FT = 2.95

    # TODO add elastic modulus? and P-V curve?
    #ε_modulus::FT = 20e6

    # TODO link to Plant Hydraulics
    "tree height `[m]`"
    height = 20
    # -Scholz et al. 2011 Book, Hydraulic Capacitance: Biophysics and Functional Significance,
    # can also be related to P50 see same book
    "tree capacitance `[kg m⁻³ Pa⁻¹]`"
    Ctree::FT = (79 + 8*height)*1e-6;

    "Use co-limitation"
    use_colim::Bool = true

    "Tree roughness `[m]`"
    z0m::FT= -999.
    # TODO should be changed later
    "tree roughness `[m]`"
    z0h::FT= -999. 
    "tree displacement height `[m]`"
    d::FT = -999.

    "leaf thickness (m)"
    dleaf::FT = 2e-3
    "Leaf area index"
    LAI::FT= 1.0
    "turbulent transfer coefficient `[m s``^2``]`"
    Cd::FT= 0.01
end








###############################################################################
#
# Meteo data set
# This struct passed the FT test
# More documentation required
#
###############################################################################
"""
    struct MeteoParams{FT}

A Placeholder for all Meteo data for the moment

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct MeteoParams{FT<:Number}
    # TODO Make this struct work first, then clean it up
    "?"
    S_down::FT = 0.0
    "?"
    L_down::FT = 0.0
    "Air temperature `[K]`"
    T_air::FT  = 280.0
    "Air vapor pressure `[Pa]`"
    e_air::FT  = 1500.0
    "Air surface pressure `[Pa]`"
    P_air::FT  = 101325.0
    
    # TODO use p_a?
    "Ambient CO₂ concentration `[µmol/mol]`"
    Ca::FT     =  400.0
    "Ambient CO₂ partial pressure `[Pa]`"
    p_a::FT    = 40.53
    "Atmospheric pressure `[Pa]`"
    p_atm::FT  = 101325.0

    "PAR `[µmol/m2/s]`"
    PAR::FT    = 600.0
    "Wind Speed `[m/s]`"
    U::FT      = 1e-6;
    "Measurement height `[m]`"
    zscreen::FT= 10.0
    "atmospheric Obukhov length"
    L::FT      = 1e6

    # parameter to define stability function for stable case
    "2 Webb correction tends to be unstable at night - suggest not using"
    stab_type_stable::FT = 1
    "?"
    ustar::FT = 1e-6

    # TODO remove this one?
    "Conversion factor from `[m/s]` to `[mol/m2/s]`"
    g_m_s_to_mol_m2_s::FT = -Inf

    # TODO remove in the future and use Pa for all the partial pressures
    "Conversion factor from mixing ratio `[μmol/mol]` to partial pressure `[Pa]`"
    ppm_to_Pa::FT = 0.1
    "Air resistance `[Pa]`"
    ra::FT    = 1e6
    "absorbed PAR `[µmol/m2/s]`"
    APAR::FT = 500.0

    # TODO use p_s instead? Note: there is a p_s in the leaf struct, remove this one?
    "CO₂ concentration at leaf surface `[µmol/mol]`"
    Cs::FT = 0.0
end