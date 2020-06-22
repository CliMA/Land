# Loading the Photosynthesis model:
using Land.Photosynthesis
# Defining our Field Type (we can easily switch between Double and Float precision this way)
const FT = Float32;

# Create a standard leaf with defualt parameters
leaf = Leaf{FT}();

# Create a standard meteo structure:
envir = AirLayer{FT}();

##?LeafParams

##?MeteoParams

##Here, we just have to worry abou the photosynthesis module, which we set here:
#mod_photo = C3FvCBPhoto()

# All modules here:
const photo_set = C3CLM(FT);

# Set APAR to 250 $\mu mol/m^2/s$
leaf.APAR = 250;
# Set temperature to 290K
leaf.T = 290;
# Applying the T-correction for all the rate constants
photo_temperature_dependence!(photo_set, leaf, envir);
@show leaf.Vcmax, leaf.Vcmax25;

@show leaf.Jmax, leaf.Jmax25;

# Specify Cc directly here in Pa
leaf.p_i = 35;
# update radiation dependent values first, like ETR
photo_radiation_dependence!(photo_set, leaf);
# update p_i dependent photosynthetic rates
photo_CO₂_dependence!(photo_set, leaf);
@show leaf.Ac;

@show leaf.Aj;

@show leaf.Ap;

@show leaf.An;

#?Leaf

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

