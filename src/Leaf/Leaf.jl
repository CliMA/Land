module Leaf
using ..PhysCon
using ..WaterVapor

export leaf_params, setLeafT!, BallBerry!, Medlyn!, setkx!, setLeafkl!, setra!,ψ

# Scaling functions for Photosynthesis temperature response and inhibition
function ft(tl, ha)
    FT = eltype(tl)
    Rgas = FT(8.31446261815324);
    T25 = FT(298.15) 
    exp(ha/(Rgas*(T25)) * (1-(T25)/tl));
end;

function fth(tl, hd, se, fc)
    FT = eltype(tl)
    Rgas = FT(8.31446261815324);
    fc / (1 + exp((-hd+se*tl)/(Rgas*tl)));
end;

function fth25(hd, se) 
    FT = eltype(hd)
    T25 = FT(298.15);
    Rgas = FT(8.31446261815324); 
    1 + exp( (-hd + se * T25) / (Rgas * T25 ));
end;


# Structure with all parameter temperature dependencies of a Leaf (largely based on Bonan's Photosynthesis model but exported as struct)
Base.@kwdef mutable struct leaf_params{TT<:Number}

    # broadband albedo and emissivity
    α::TT  = -999;                           # broadband shortwave albedo - absurd value to make sure we initialize correctly
    ε::TT  = -999;                           # longwave emissivity

    # thermal characteristics
    LMA::TT         = 100e-3;                # DRY leaf mass area, NOT MEAN VALUE WHICH INCLUDES RELATIVE WATER CONTENT - kg/m2 - density times thickness
    RWC::TT         = 0.8;                   # leaf relative water content
    Cleaf::TT       = 1000.0;                # leaf specific heat (J/kg/K)

    # pore traits
    #pore_density::TT   = 100.0  / (1e-6);     # pore density in pores/m^2
    Chloroplast_rel_volume::TT  = 2.0/100.0;        # hloroplast volume airspace (m^3 per pore) - http://kirschner.med.harvard.edu/files/bionumbers/A%20semi%20quantitative%20analysis%20of%20the%20photosynthetic%20system%20in%20an%20'average'%20C3%20plant%20leaf.pdf

    # Rate constants (arbitrary units here, only relative rates are important):
    Kf::TT = 0.05;                           # Rate constant for fluorescence (might try to fit this eventually)
    Kd::TT = 0.85;                           # Rate constant for thermal dissipation
    Kp::TT = 4.0;                            # Rate constant for photochemistry (all reaction centers open)
    Kn::TT = 0;                              # NPQ rate constant (initially zero)
    Kn_ss::TT = 0;                           # Kn in steady state
    maxPSII::TT = Kp/(Kp+Kf+Kd);             # Max PSII yield (Kn=0, all RC open)
    "Flexas et al derived Fluorescence model params"
    Knparams = [5.01, 1.93, 10.0]
    Fo::TT = 0;     # dark-adapted fluorescence yield Fo,0
    Fo′::TT = 0;     # light-adapted fluorescence yield in the dark Fo
    Fm::TT = 0;     # light-adapted fluorescence yield Fm
    Fm′::TT = 0;    # dark-adapted fluorescence yield Fm
    ϕs::TT = 0;     # steady-state (light-adapted) yield Ft (aka Fs)
    eta::TT = 0;
    qQ::TT = 0;     # photochemical quenching
    qE::TT = 0;     # non-photochemical quenching
    NPQ::TT = 0;
    "Dynamic Steady state boolean"
    dynamic_state::Bool = false


    effcon::TT = 1/5                       # [mol CO2/mol e-]  number of CO2 per electrons - typically 1/5 for C3 and 1/6 for C4, i.e. about 1/10 for both PSI and PSII!
    CO2_per_electron::TT = 1/5             # Similar to above but accounting for Photorespiration as well (which "steals" electrons)
    Ja::TT = 0                               # Actual electron transport rate (needed to compute fluorescence later on)

    C3::Bool = true                           # C3 (or C4) plant
    use_colim::Bool = true                    # Use co-limitation
    gstyp::Int = 1                            # stomatal conductance type (1=Ball Berry, 0 = Medlyn)
    o₂::TT = 0.209e3;                         # Standard O2 in mmol/mol
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

    gleaf::TT = 0.1;                           # total leaf conductance (μmol/m2/s)
    gm::TT = Inf;                              # Mesophyll conductance (μmol/m2/s/)
    gs::TT = 0.1;                              # Stomatal conductance (μmol/m2/s); just using some prior
    gs_ss::TT = 0.1;                           # Steady state Stomatal conductance (μmol/m2/s);

    #---------------------------------------------------------------------
    # Kc, Ko, cp at 25C: Bernacchi et al (2001) Plant, Cell Environment 24:253-259
    # Derive sco from cp with o2=0.209 mol/mol and re-calculate Γ to allow
    # variation in O2
    #---------------------------------------------------------------------
    Kc_25::TT = 404.9;                         # Michaelis-Menten constant for CO2 at 25C - μmol/mol
    Ko_25::TT = 278.4;                         # Michaelis-Menten constant for O2 at 25C  - mmol/mol
    Γ_25_::TT = 42.75;                         # CO2 compensation point at 25C μmol/mol
    sco::TT = 0.5 * 0.209 / (Γ_25_ * 1.e-06);  # Γ_25 (μmol/mol) -> (mol/mol)
    Γ_25::TT = 0.5 * o₂ / sco * 1000.;         # O2 is mmol/mol. Multiply by 1000 for μmol/mol

    #---------------------------------------------------------------------
    # Activation energy:
    # Bernacchi et al (2001) Plant, Cell Environment 24:253-259
    # Bernacchi et al (2003) Plant, Cell Environment 26:1419-1430
    # Acclimation from: Kattge and Knorr (2007) Plant, Cell Environment 30:1176-1190
    #---------------------------------------------------------------------
    Kcha::TT    = 79430.;                      # Activation energy for Kc (J/mol)
    Koha::TT    = 36380.;                      # Activation energy for Ko (J/mol)
    cpha::TT    = 37830.;                      # Activation energy for cp (J/mol)
    Vcmaxha::TT = 65330.;                      # Activation energy for Vcmax (J/mol)
    Jmaxha::TT  = 43540.;                      # Activation energy for Jmax (J/mol)
    rdha::TT    = 46390.;                      # Activation energy for Rd (J/mol)

    #---------------------------------------------------------------------
    # High temperature deactivation:
    # Leuning (2002) Plant, Cell Environment 25:1205-1210
    # The factor "c" scales the deactivation to a value of 1.0 at 25C
    # Acclimation from: Kattge and Knorr (2007) Plant, Cell Environment 30:1176-1190
    #---------------------------------------------------------------------
    Vcmaxhd::TT = 150000.;                    # Deactivation energy for Vcmax (J/mol)
    Jmaxhd::TT  = 150000.;                    # Deactivation energy for Jmax (J/mol)
    rdhd::TT    = 150000.;                    # Deactivation energy for Rd (J/mol)

    Vcmaxse::TT = 490.;                       # Entropy term for Vcmax (J/mol/K)
    Jmaxse::TT  = 490.;                       # Entropy term for Jmax (J/mol/K)
    rdse::TT    = 490.;                       # Entropy term for Rd (J/mol/K)

    Vcmaxc::TT = fth25(Vcmaxhd, Vcmaxse);    # Scaling factor for high temperature inhibition (25 C = 1.0)
    Jmaxc::TT  = fth25(Jmaxhd, Jmaxse);      # Scaling factor for high temperature inhibition (25 C = 1.0)
    rdc::TT    = fth25(rdhd, rdse);          # Scaling factor for high temperature inhibition (25 C = 1.0)

    # Initialize at standard T:
    Kc::TT  = Kc_25;
    Ko::TT  = Ko_25;
    Γstar::TT  = Γ_25;

    # Use some standard values first
    Vcmax25::TT = 80.0;                        # Leaf maximum carboxylation rate at 25C for canopy layer (μmol/m2/s)
    Jmax25::TT  = 1.97*Vcmax25;                # C3 - maximum electron transport rate at 25C for canopy layer (μmol/m2/s)
    Rd25::TT    = 0.7605;                      # Leaf respiration rate at 25C for canopy layer (μmol CO2/m2/s)
    Vcmax::TT  = Vcmax25;
    Jmax::TT   = Jmax25;
    Rdleaf::TT = Rd25;


    Cc::TT = 400.0;                               # CO2 concentration in chloroplast [ppm]
    RH::TT = 100.0;                               # relative humidity at the surface of the leaf
    VPD::TT = 0.1;                                # VPD gradient ACROSS THE LEAF interface - NOT the weather station one 8-)

    # to be computed
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
    l.Kc     = l.Kc_25          * ft(l.T, l.Kcha);
    l.Ko     = l.Ko_25          * ft(l.T, l.Koha);
    l.Γstar  = l.Γ_25           * ft(l.T, l.cpha);
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

include("leaf_photosynthesis.jl")
include("math_tools.jl")
include("leaf_energy_water_balance.jl")

end #Module
