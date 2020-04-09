
export leaf_params, setLeafT!, BallBerry!, Medlyn!, setkx!, setra!,ψ

# Scaling functions for Photosynthesis temperature response and inhibition
ft(tl, ha) = exp(ha/(physcon.Rgas*(physcon.tfrz+25)) * (1-(physcon.tfrz+25)/tl));
fth(tl, hd, se, fc) = fc / (1 + exp((-hd+se*tl)/(physcon.Rgas*tl)));
fth25(hd, se) = 1.0 + exp( (-hd + se * (physcon.tfrz+25.)) / (physcon.Rgas * (physcon.tfrz+25.)) );

# Structure with all parameter temperature dependencies of a Leaf (largely based on Bonan's Photosynthesis model but exported as struct)
@with_kw mutable struct leaf_params{TT<:Number}

    # broadband albedo and emissivity
    α::TT  = -999;                           # broadband shortwave albedo - absurd value to make sure we initialize correctly
    ε::TT  = -999;                           # longwave emissivity

    # thermal characteristics
    LMA::TT         = -999;                  # leaf mass area - kg/m2 - density times thickness
    c_leaf::TT      = -999;                  # leaf heat capacity - J/kg

    # turbulence
    ra::TT          = 20;                    # leaf aerodynamic resistance (s/m)

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
    g0::TT = 0.01;                          # Ball-Berry minimum leaf conductance, unstressed (mol H2O/m2/s)
    g1::TT = 9;                             # Ball-Berry slope of conductance-photosynthesis relationship, unstressed
    vpd_min = 100.0                            # Min VPD for Medlyn model (blows up otherwise?)
    # Add already here: Mesophyll conductance (Inf for now):
    gm::TT = Inf;                              # Mesophyll conductance (μmol/m2/s/Pa)
    gs::TT = 0.1;                              # Stomatal conductance (μmol/m2/s/Pa); just using some prior
    gs_ss::TT = 0.1;                           # Steady state Stomatal conductance (μmol/m2/s/Pa);

    #---------------------------------------------------------------------
    # kc, ko, cp at 25C: Bernacchi et al (2001) Plant, Cell Environment 24:253-259
    # Derive sco from cp with o2=0.209 mol/mol and re-calculate Γ to allow
    # variation in O2
    #---------------------------------------------------------------------
    kc_25::TT = 404.9;                         # Michaelis-Menten constant for CO2 at 25C μmol/mol
    ko_25::TT = 278.4;                         # Michaelis-Menten constant for O2 at 25C  mmol/mol
    Γ_25_::TT = 42.75;                         # CO2 compensation point at 25C μmol/mol
    sco::TT = 0.5 * 0.209 / (Γ_25_ * 1.e-06);  # Γ_25 (μmol/mol) -> (mol/mol)
    Γ_25::TT = 0.5 * o₂ / sco * 1000.;         # O2 is mmol/mol. Multiply by 1000 for μmol/mol

    #---------------------------------------------------------------------
    # Activation energy:
    # Bernacchi et al (2001) Plant, Cell Environment 24:253-259
    # Bernacchi et al (2003) Plant, Cell Environment 26:1419-1430
    # Acclimation from: Kattge and Knorr (2007) Plant, Cell Environment 30:1176-1190
    #---------------------------------------------------------------------
    kcha::TT    = 79430.;                      # Activation energy for kc (J/mol)
    koha::TT    = 36380.;                      # Activation energy for ko (J/mol)
    cpha::TT    = 37830.;                      # Activation energy for cp (J/mol)
    vcmaxha::TT = 65330.;                      # Activation energy for Vcmax (J/mol)
    jmaxha::TT  = 43540.;                      # Activation energy for Jmax (J/mol)
    rdha::TT    = 46390.;                      # Activation energy for Rd (J/mol)

    #---------------------------------------------------------------------
    # High temperature deactivation:
    # Leuning (2002) Plant, Cell Environment 25:1205-1210
    # The factor "c" scales the deactivation to a value of 1.0 at 25C
    # Acclimation from: Kattge and Knorr (2007) Plant, Cell Environment 30:1176-1190
    #---------------------------------------------------------------------
    vcmaxhd::TT = 150000.;                    # Deactivation energy for Vcmax (J/mol)
    jmaxhd::TT  = 150000.;                    # Deactivation energy for Jmax (J/mol)
    rdhd::TT    = 150000.;                    # Deactivation energy for Rd (J/mol)

    vcmaxse::TT = 490.;                       # Entropy term for Vcmax (J/mol/K)
    jmaxse::TT  = 490.;                       # Entropy term for Jmax (J/mol/K)
    rdse::TT    = 490.;                       # Entropy term for Rd (J/mol/K)

    vcmaxc::TT = fth25(vcmaxhd, vcmaxse);    # Scaling factor for high temperature inhibition (25 C = 1.0)
    jmaxc::TT  = fth25(jmaxhd, jmaxse);      # Scaling factor for high temperature inhibition (25 C = 1.0)
    rdc::TT    = fth25(rdhd, rdse);          # Scaling factor for high temperature inhibition (25 C = 1.0)

    # Initialize at standard T:
    kc::TT  = kc_25;
    ko::TT  = ko_25;
    Γstar::TT  = Γ_25;

    # Use some standard values first
    vcmax25::TT = 80.0;                        # Leaf maximum carboxylation rate at 25C for canopy layer (μmol/m2/s)
    jmax25::TT  = 1.97*vcmax25;                # C3 - maximum electron transport rate at 25C for canopy layer (μmol/m2/s)
    rd25::TT    = 0.7605;                      # Leaf respiration rate at 25C for canopy layer (μmol CO2/m2/s)
    vcmax::TT  = vcmax25;
    jmax::TT   = jmax25;
    rdleaf::TT = rd25;

    #
    Ci::TT = 0.0;                               # CO2 concentration in mesophyll/internal

    # to be computed
    je::TT = 0.0 ;                              # electron transport rate





    # tree/leaf traits
    height      = 20.;                            # tree height (m)
    z0m         = 0.1*height;                     # tree roughness (m)
    d           = 2/3*height;                     # tree displacement height (m)

    dleaf       = 2e-3;                           # leaf thickness (m)
    Cd          = 0.01;                           # m/sqrt(s) turbulent transfer coefficient

    # plant hydraulics
    psi_l::TT   = -1.5e6;                         # leaf water potential (Pa)

    psi_l50::TT = -1.75e6;                        # leaf water potential at 50% drop in conductivity (Pa)
    kmax::TT    = 4e-8                            # maximum leaf-xylem conductivity (m/s)
    kx::TT      = kmax                            # actual xylem conductivity (m/s)
    ck::TT      = 2.95;                           # slope of Weibull curve
    #ε_modulus::TT = 20e6;                        # elastic modulus - for later
    Ctree::TT   = (79 + 8*height)*1e-6;           # tree capacitance (kg m−3 Pa−1)
                                                  # -Scholz et al. 2011 Book, Hydraulic Capacitance: Biophysics and Functional Significance,
                                                  # can also be related to P50 see same book
end

# Set Leaf rates with vcmax, jmax and rd at 25C as well as actual T here:
# For some reason, this is slow and allocates a lot, can be improved!!
"Set Leaf rates with vcmax, jmax and rd at 25C as well as actual T here"
function setLeafT!(l::leaf_params,  T)
    l.T      = T;
    l.kc     = l.kc_25          * ft(T, l.kcha);
    l.ko     = l.ko_25          * ft(T, l.koha);
    l.Γstar  = l.Γ_25           * ft(T, l.cpha);
    l.vcmax   = l.vcmax25 * ft(T, l.vcmaxha) * fth(T, l.vcmaxhd, l.vcmaxse, l.vcmaxc);
    l.jmax    = l.jmax25  * ft(T, l.jmaxha)  * fth(T, l.jmaxhd, l.jmaxse, l.jmaxc);
    l.rdleaf  = l.rd25    * ft(T, l.rdha)    * fth(T, l.rdhd, l.rdse, l.rdc);
    (l.esat, l.desat) .= SatVap(T);
    # l.kd = max(0.8738,  0.0301*(T-273.15)+ 0.0773); # Can implement that later.
end

# Ball-Berry stomatal conductance model:
function BallBerry!(Cs, RH, A, l::leaf_params)
  #  Cs  : CO2 at leaf surface
  #  RH  : relative humidity
  #  A   : Net assimilation in 'same units of CO2 as Cs'/m2/s
  #  minCi : minimum Ci as a fraction of Cs (in case RH is very low?)
  #  Ci_input : will calculate gs if A is specified.

  l.gs = l.g1 * A * RH/Cs  + l.g0;
end # function

# Medlyn stomatal conductance model:
function Medlyn!(Cs, VPD, A, l::leaf_params)
  #  Cs  : CO2 at leaf surface
  #  VPD  : vapor pressure deficit - Medlyn model
  #  A   : Net assimilation in 'same units of CO2 as Cs'/m2/s
  #  minCi : minimum Ci as a fraction of Cs (in case RH is very low?)
  #  Ci_input : will calculate gs if A is specified.

  #print(A,Cs,l.g0,l.g1,l.gs)
  l.gs = (1+l.g1/sqrt(VPD)) * A /Cs  + l.g0;
  #print(l.gs)
end # function

function setkx!(l::leaf_params, psis, psi_l) # set hydraulic conductivity
    l.kx = l.kmax * IntWeibull(psis,psi_l,l.psi_l50,l.ck); # kmax . int_psis^psil k(x)dx = kmax . IntWeibull(psil);
end
