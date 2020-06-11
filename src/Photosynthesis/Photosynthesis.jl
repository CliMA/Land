module Photosynthesis

using DocStringExtensions
using Parameters

using ..LandParameters

@unpack GAS_R, K_25 = LandParameters

export AbstractPhotoModelParaSet,
       C3ParaSet,
       C4ParaSet,
       FluoParaSet,
       C3Bernacchi,
       C3CLM,
       C4CLM,
       FluorescenceFlexas,
       arrhenius_correction,
       an_ag_r_from_pi,
       an_ag_r_pi_from_gsc,
       get_Γ_star




include("types.jl"     )
include("parasets.jl"  )
include("photomodel.jl")




###############################################################################
#
# Functions in development
# Test and document the functions before merging them to this file
#
###############################################################################
include("Photosynthesis_in_development.jl")




end