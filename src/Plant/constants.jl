# set the default floating type, Float64 for now. Should be the same as in the RT module
const FT = Float64

# load the constant gravity value
struct EarthParameterSet <: AbstractEarthParameterSet end
const earth   = EarthParameterSet()
const gravity = FT( CLIMAParameters.Planet.grav(earth) )
const GAS_R   = FT( gas_constant() )

# temperature constants
const K_0    = FT(273.15)
const K_25   = FT(298.15)




# more pending, may need to migrate this to CLIMAParameters in the future
const ρ_H₂O  = FT(998.0 )    # Kg m⁻³ | @ 293.15 K, may need to calculate it from temperature in the future
const svp_a0 = FT(611.0 )    # Pa     | Used for satrurated vapor pressure
const svp_e0 = FT(17.502)    #        | Used for satrurated vapor pressure
const svp_e1 = FT(240.97)    #        | Used for satrurated vapor pressure

const slh_a0 = FT(2500.8)    # | Used for specific latent heat
const slh_a1 = FT(-2.36 )    # | Used for specific latent heat
const slh_a2 = FT(1.6e-3)    # | Used for specific latent heat
const slh_a3 = FT(-6e-5 )    # | Used for specific latent heat

const vis_a = FT(1.856E-14)    # Pa s | Used for viscosity
const vis_b = FT(4209     )    # K    | Used for viscosity
const vis_c = FT(0.04527  )    # K⁻¹  | Used for viscosity
const vis_d = FT(-3.376E-5)    # K⁻²  | Used for viscosity

const st_k  = FT(2.1E-7)    # J K⁻¹ mol⁻(2/3) | Used for surface tension
const st_v  = FT(18.0  )    # ml mol⁻¹        | Used for surface tension
const st_tc = FT(647.0 )    # K               | Used for surface tension




# constants for photosynthesis model, may need to merge with the Photosynthesis in the RT module
# constants for Vcmax
const PS_V_HA = FT(73637.0)
const PS_V_HD = FT(149252.0)
const PS_V_SV = FT(486.0)
const PS_V_C  = 1 + exp( (PS_V_SV*K_25 - PS_V_HD) / (GAS_R*K_25) )

# constants for Jmax and J
const PS_J_HA = FT(50300.0)
const PS_J_HD = FT(152044.0)
const PS_J_SV = FT(495.0)
const PS_J_C  = 1 + exp( (PS_J_SV*K_25 - PS_J_HD) / (GAS_R*K_25) )

const PS_J_QY = FT(0.3)    # quantum yield
const PS_J_CR = FT(0.9)    # curvature factor

const PS_RV_RATIO = FT(0.015 )
const PS_R_Q10_E  = FT(2.0   )
const PS_R_EXP_K  = FT(1.3   )
const PS_R_EXP_T  = FT(328.15)

const PS_KC_25  = FT(41.01637)
const PS_KC_Q10 = FT(2.1     )
const PS_KO_25  = FT(28201.92)
const PS_KO_Q10 = FT(1.2     )
