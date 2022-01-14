###############################################################################
#
# Fluorescence model parameter set
# Data source: van der Tol et al. (2014)
# Models of fluorescence and photosynthesis for interpreting measurements of
#     solar-induced chlorophyll fluorescence
#
###############################################################################
""" [`FluoParaSet`](@ref) type parameter set using all data """
FluorescenceVanDerTol(FT) = FluoParaSet{FT}(2.48, 2.83, 0.114);

""" [`FluoParaSet`](@ref) type parameter set using Flexas's data (drought) """
FluorescenceVanDerTolDrought(FT) = FluoParaSet{FT}(5.01, 1.93, 10.0);
