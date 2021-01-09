#=
###############################################################################
#
# Leaf photosynthesis-related parameter set
# This struct passed the FT test
# This struct is documented in the Leaf page
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
    "Kn in steady state"
    Kn_ss::FT = 0

    # Fluorescence-related
    "Dynamic Steady state `[Bool]`"
    dynamic_state::Bool = false

    # number of CO₂ per electrons - typically 1/5 for C3 and 1/6 for C4, i.e. about 1/10 for both PSI and PSII!
    "Efficiency of photosynthesis light reactions`[mol CO₂/mol e-]`"
    effcon::FT = 1/5

    "Derivative of saturation vapor pressure at leaf temperature wrt T `[Pa/K]`"
    desat::FT = 0.0

    "Steady state Stomatal conductance `[mol H₂O m⁻² s⁻¹]`"
    gs_ss::FT = 0.1

    # TODO Update this with Plant Hydraulics
    "Sap flow `[?]`"
    Sap::FT = 0;
    "relative humidity at the surface of the leaf `[-]`"
    RH::FT = 100.0

    # VPD gradient ACROSS THE LEAF interface - NOT the weather station one 8-)
    "Vapor pressure difference across leaf interface `[Pa]`"
    VPD::FT = 0.1

    # TODO add elastic modulus? and P-V curve?
    #ε_modulus::FT = 20e6

    # TODO link to Plant Hydraulics
    "tree height `[m]`"
    height = 20
    # -Scholz et al. 2011 Book, Hydraulic Capacitance: Biophysics and Functional Significance,
    # can also be related to P50 see same book
    "tree capacitance `[kg m⁻³ Pa⁻¹]`"
    Ctree::FT = (79 + 8*height)*1e-6;

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

    "Measurement height `[m]`"
    zscreen::FT= 10.0
    "atmospheric Obukhov length"
    L::FT      = 1e6

    # parameter to define stability function for stable case
    "2 Webb correction tends to be unstable at night - suggest not using"
    stab_type_stable::FT = 1
    "?"
    ustar::FT = 1e-6

    "Air resistance `[Pa]`"
    ra::FT    = 1e6
end
=#
