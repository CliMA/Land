using CLIMAParameters

# set the default floating type, Float32 for now. Should be the same as in the RT module
FT = Float64

# more pending, may need to migrate this to CLIMAParameters in the future
ρ_H₂O = FT(998.0) # Kg m⁻³ | @ 293.15 K, may need to calculate it from temperature in the future

# load the constant gravity value
struct EarthParameterSet <: AbstractEarthParameterSet end
earth   = EarthParameterSet()
const gravity = CLIMAParameters.Planet.grav(earth)
