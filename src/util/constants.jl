#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jul-08: add struct for universal constants
#     2022-Jul-08: add field and thus wrapper for thermal conductivity of liquid water
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save universal constants.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct UniversalConstants
    "Avogadro's number `[molecule mol⁻¹]`"
    AVOGADRO::Float64 = 6.02214076e23
    "Isobaric specific heat of liquid water `[J kg⁻¹ K⁻¹]`"
    CP_L::Float64 = 4181
    "Isobaric specific heat of water vapor `[J kg⁻¹ K⁻¹]`"
    CP_V::Float64 = 1859
    "Universal gas constant `[J mol⁻¹ K⁻¹]`"
    GAS_R::Float64 = 8.3144598
    "Gravity of the Earth `[m s⁻²]`"
    GRAVITY::Float64 = 9.81
    "Planck constant `[m² kg s⁻¹]`"
    H_PLANCK::Float64 = 6.626e-34
    "Boltzmann constant `[m² kg s⁻² K⁻¹]`"
    K_BOLTZMANN::Float64 = 1.381e-23
    "Stefan-Boltzmann constant `[W m⁻² K⁻⁴]`"
    K_STEFAN::Float64 = 5.670e-8
    "Von Karman constant `[-]`"
    K_VON_KARMAN::Float64 = 0.4
    "Latent heat vaporization at ``T_0`` `[K kg⁻¹]`"
    LH_V0::Float64 = 2.5008e6
    "Light speed in vacuum `[m s⁻¹]`"
    LIGHT_SPEED::Float64 = 2.99792458e8
    "Molar mass of dry air `[kg mol⁻¹]`"
    M_DRYAIR::Float64 = 28.97e-3
    "Molar mass of water `[kg mol⁻¹]`"
    M_H₂O::Float64 = 18.01528e-3
    "Mean atmospheric pressure at sea level `[Pa]`"
    P_ATM::Float64 = 1.01325e5
    "Water vapor pressure at triple temperature `[Pa]`"
    PRESS_TRIPLE::Float64 = 611.657
    "Freezing temperature of water `[K]`"
    T_0::Float64 = 273.15
    "Triple temperature of water `[K]`"
    T_TRIPLE::Float64 = 273.16
    "Mean number of days per year [day]"
    YEAR_D::Float64 = 365.2422222
    "Thermal conductivity of water `[W m⁻¹ K⁻¹]`"
    Λ_THERMAL_H₂O::Float64 = 0.57
    "Density of liquid water `[kg m⁻³]`"
    ρ_H₂O::Float64 = 1000
end


#######################################################################################################################################################################################################
#
# Changes to these wrappers
# General
#     2022-Jul-08: add wrapper functions
#
#######################################################################################################################################################################################################
const UNIVERSAL_CONSTANTS = UniversalConstants();

""" Avogadro's number `[molecule mol⁻¹]` """
AVOGADRO(FT=Float64) = FT(UNIVERSAL_CONSTANTS.AVOGADRO);

""" Isobaric specific heat of dry air `[J kg⁻¹ K⁻¹]` """
CP_D(FT=Float64) = GAS_R(FT) / M_DRYAIR(FT) / FT(3.5);

""" Isobaric specific heat of dry air `[J mol⁻¹ K⁻¹]` """
CP_D_MOL(FT=Float64) = GAS_R(FT) / FT(3.5);

""" Isobaric specific heat of liquid water `[J kg⁻¹ K⁻¹]` """
CP_L(FT=Float64) = FT(UNIVERSAL_CONSTANTS.CP_L);

""" Isobaric specific heat of water vapor `[J kg⁻¹ K⁻¹]` """
CP_V(FT=Float64) = FT(UNIVERSAL_CONSTANTS.CP_V);

""" Universal gas constant `[J mol⁻¹ K⁻¹]` """
GAS_R(FT=Float64) = FT(UNIVERSAL_CONSTANTS.GAS_R);

""" Gravity of the Earth `[m s⁻²]` """
GRAVITY(FT=Float64) = FT(UNIVERSAL_CONSTANTS.GRAVITY);

""" Planck constant `[m² kg s⁻¹]` """
H_PLANCK(FT=Float64) = FT(UNIVERSAL_CONSTANTS.H_PLANCK);

""" Boltzmann constant `[m² kg s⁻² K⁻¹]` """
K_BOLTZMANN(FT=Float64) = FT(UNIVERSAL_CONSTANTS.K_BOLTZMANN);

""" Stefan-Boltzmann constant `[W m⁻² K⁻⁴]` """
K_STEFAN(FT=Float64) = FT(UNIVERSAL_CONSTANTS.K_STEFAN);

""" Von Karman constant `[-]` """
K_VON_KARMAN(FT=Float64) = FT(UNIVERSAL_CONSTANTS.K_VON_KARMAN);

""" Latent heat vaporization at ``T_0`` `[K kg⁻¹]` """
LH_V0(FT=Float64) = FT(UNIVERSAL_CONSTANTS.LH_V0);

""" Light speed in vacuum `[m s⁻¹]` """
LIGHT_SPEED(FT=Float64) = FT(UNIVERSAL_CONSTANTS.LIGHT_SPEED);

""" Molar mass of dry air `[kg mol⁻¹]` """
M_DRYAIR(FT=Float64) = FT(UNIVERSAL_CONSTANTS.M_DRYAIR);

""" Molar mass of water `[kg mol⁻¹]` """
M_H₂O(FT=Float64) = FT(UNIVERSAL_CONSTANTS.M_H₂O);

""" Mean atmospheric pressure at sea level `[Pa]` """
P_ATM(FT=Float64) = FT(UNIVERSAL_CONSTANTS.P_ATM);

""" Water vapor pressure at triple temperature `[Pa]` """
PRESS_TRIPLE(FT=Float64) = FT(UNIVERSAL_CONSTANTS.PRESS_TRIPLE);

""" Gas constant water vapor `[J kg⁻¹ K⁻¹]` """
R_V(FT=Float64) = GAS_R(FT) / M_H₂O(FT);

""" Gas constant times 298.15 K `[J mol⁻¹]` """
RT_25(FT=Float64) = GAS_R(FT) * T_25(FT);

""" Freezing temperature of water `[K]` """
T_0(FT=Float64) = FT(UNIVERSAL_CONSTANTS.T_0);

""" Kelvin temperature at 25 Celcius `[K]` """
T_25(FT=Float64) = T_0(FT) + 25;

""" Triple temperature of water `[K]` """
T_TRIPLE(FT=Float64) = FT(UNIVERSAL_CONSTANTS.T_TRIPLE);

""" Molar volume of liqiud water """
V_H₂O(FT=Float64) = M_H₂O(FT) / ρ_H₂O(FT);

""" Mean number of days per year [day] """
YEAR_D(FT=Float64) = FT(UNIVERSAL_CONSTANTS.YEAR_D);

""" Thermal conductivity of water `[W m⁻¹ K⁻¹]` """
Λ_THERMAL_H₂O(FT=Float64) = FT(UNIVERSAL_CONSTANTS.Λ_THERMAL_H₂O);

""" Density of liquid water `[kg m⁻³]` """
ρ_H₂O(FT=Float64) = FT(UNIVERSAL_CONSTANTS.ρ_H₂O);

""" Density of water times gravity `[MPa m⁻¹]` """
ρg_MPa(FT=Float64) = ρ_H₂O(FT) * GRAVITY(FT) * FT(1e-6);
