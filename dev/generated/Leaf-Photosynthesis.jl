# Loading the Photosynthesis model:
using Revise
using Land.Leaf
using Land.Photosynthesis
# Defining our Field Type (we can easily switch between Double and Float precision this way)
const FT = Float32;

# Create a standard leaf with defualt parameters
leaf = LeafParams{FT}();

# Create a standard meteo structure:
met = MeteoParams{FT}();

#?LeafParams

#?MeteoParams

##Here, we just have to worry abou the photosynthesis module, which we set here:
#mod_photo = C3FvCBPhoto()

# All modules here:
mods = C3CLM(FT)

# Set APAR to 250 $\mu mol/m^2/s$
leaf.APAR = 250;
# Set temperature to 290K
leaf.T = 290;
# Applying the T-correction for all rate constants:
Leaf.update_leaf_TD!(mods, leaf)
@show leaf.Vcmax
@show leaf.Jmax
# Specify Cc directly here in ppm (will be converted to Pa internally)
leaf.Cc = 350;

Aj = Leaf.light_limited_rate!(mods, leaf, leaf.APAR)

Ac = Leaf.rubisco_limited_rate!(mods, leaf)

Ap = Leaf.product_limited_rate!(mods, leaf)

#?C3FvCBPhoto

#?C4CollatzPhoto

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

