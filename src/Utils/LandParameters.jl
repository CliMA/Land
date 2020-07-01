module LandParameters

using CLIMAParameters
Planet       = CLIMAParameters.Planet
SubgridScale = CLIMAParameters.SubgridScale




###############################################################################
#
# Create an Earth planet
# This struct is documented in the Utils page
#
###############################################################################
""" Define a local struct inherited from AbstractEarthParameterSet  """
struct EarthParameterSet <: AbstractEarthParameterSet end
""" Earth as EarthParameterSet """
const EARTH = EarthParameterSet()








###############################################################################
#
# GLOBAL constants from CLIMAParameters
# These variables are documented in the Utils page
#
###############################################################################
""" Avogadro's number"""
const AVOGADRO         = avogad()
""" Isobaric specific heat of dry air `[J kg⁻¹ K⁻¹]` """
const CP_D             = Planet.cp_d(EARTH)
""" Isobaric specific heat of ice `[J kg⁻¹ K⁻¹]` """
const CP_I             = Planet.cp_i(EARTH)
""" Isobaric specific heat of liquid water `[J kg⁻¹ K⁻¹]` """
const CP_L             = Planet.cp_l(EARTH)
""" Isobaric specific heat of water vapor vapor `[J kg⁻¹ K⁻¹]` """
const CP_V             = Planet.cp_v(EARTH)
""" Universal gas constant `[J K⁻¹ mol⁻¹]` """
const GAS_R            = gas_constant()
""" Gravitational acceleration `[m s⁻²]` """
const GRAVITY          = Planet.grav(EARTH)
""" Planck constant `[m² kg s⁻¹]` """
const H_PLANCK         = h_Planck()
""" Temperature at 0 Celcius `[K]` """
const K_0              = Planet.T_freeze(EARTH)
""" Temperature at 25 Celcius `[K]` """
const K_25             = Planet.T_freeze(EARTH) + 25
""" Stefan-Boltzmann constant (W m⁻² K⁻⁴) """
const K_BOLTZMANN      = k_Boltzmann()
""" Latent heat sublimation at ``T_{triple}`` `[J kg⁻¹]` """
const LH_S0            = Planet.LH_s0(EARTH)
""" Latent heat vaporization at ``T_{triple}`` `[J kg⁻¹]` """
const LH_V0            = Planet.LH_v0(EARTH)
""" Light speed `[m s⁻¹]` """
const LIGHT_SPEED      = light_speed()
""" Mole mass of dry air `[kg mol⁻¹]` """
const MOLMASS_DRYAIR   = Planet.molmass_dryair(EARTH)
""" Mole mass of water `[kg mol⁻¹]` """
const MOLMASS_WATER    = Planet.molmass_water(EARTH)
""" Mean sea level pressure `[Pa]` """
const P_ATM            = Planet.MSLP(EARTH)
""" Triple point vapor pressure `[Pa]` """
const PRESS_TRIPLE     = Planet.press_triple(EARTH)
""" Ideal gas constant for dry air [J kg⁻¹ K⁻¹] """
const R_D              = Planet.R_d(EARTH)
""" Ideal gas constant for water vapor [J kg⁻¹ K⁻¹] """
const R_V              = Planet.R_v(EARTH)
""" Product of GAS_R and K_25 """
const RK_25            = GAS_R * K_25
""" Triple point temperature `[K]` """
const T_TRIPLE         = Planet.T_triple(EARTH)
""" von Karman constant """
const VON_KARMAN_CONST = SubgridScale.von_karman_const(EARTH)
""" Ratio of molecular mass of water and dry air """
const WATER_AIR_MRATIO = 1 / Planet.molmass_ratio(EARTH)
""" Days per year """
const YEAR_D           = 365.2422222
""" Density of water `[kg m⁻³]`"""
const ρ_H₂O            = Planet.ρ_cloud_liq(EARTH)




end # module
