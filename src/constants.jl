# define the constants for Earth.jl
struct EarthParameterSet <: AbstractEarthParameterSet end
const EARTH         = EarthParameterSet()
const GAS_R         = gas_constant()
const GRAVITY       = Planet.grav(EARTH)
const K_0           = Planet.T_freeze(EARTH)
const K_25          = K_0 + 25
const K_BOLTZMANN   = k_Boltzmann()
const MOLMASS_WATER = Planet.molmass_water(EARTH)
const P_ATM         = Planet.MSLP(EARTH)
const RK_25         = GAS_R * K_25
const YEAR_D        = 365.2422222


const KG_H_2_MOL_S  = 1000 / MOLMASS_WATER / 3600


