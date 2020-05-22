# Loading the Photosynthesis model:
using Land.Photosynthesis
# Defining our Field Type (we can easily switch between Double and Float precision this way)
const FT = Float32;

# Create a standard leaf with defualt parameters
leaf = leaf_params{FT}();

# Create a standard meteo structure:
met = meteo{FT}();

?leaf_params

?meteo

##Here, we just have to worry abou the photosynthesis module, which we set here:
mod_photo = C3FvCBPhoto()

# Set APAR to 250 $\mu mol/m^2/s$
leaf.APAR = 250;
# Set temperature to 290K
leaf.T = 290;
# Applying the T-correction for all rate constants:
Photosynthesis.set_leaf_temperature!(mods, leaf)
@show leaf.Vcmax
@show leaf.Jmax
# Specify Cc directly here in ppm (will be converted to Pa internally)
leaf.Cc = 350;

Aj = Photosynthesis.light_limited_rate!(mod_photo, leaf, met, leaf.APAR)

Ac = Photosynthesis.rubisco_limited_rate!(mod_photo, leaf, met)

Ap = Photosynthesis.product_limited_rate!(mod_photo, leaf, met)

?C3FvCBPhoto

?C4CollatzPhoto

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

