module Photosynthesis

using DocStringExtensions
using Parameters

using ..PhysCon
using ..WaterVapor

export leaf_params, setLeafT!, BallBerry!, Medlyn!, Gentine!, setkx!, setLeafkl!, setra!,ψ

export fluxes, meteo, ψ_h, ψ_m, setra!, setRoughness!, LeafPhotosynthesis!

export AbstractLeafRespiration, RespirationCLM, RespirationBernacchi

export AbstractPhotosynthesisModel, C3FvCBPhoto, C3FvCBPhotoGs, C3FvCBPhotoATP, C4CollatzPhoto

export AbstractJmax, AbstractVmax, AbstractVcmax, AbstractVpmax

export AbstractMM 

export JmaxCLM, VcmaxCLM, JmaxBernacchi, VcmaxBernacchi, Vpmax, MM_CLM 

export AbstractFluorescenceModel, FlexasTolBerry_fluorescence

export rubisco_limited_rate, light_limited_rate

include("leaf.jl")
include("leaf_photosynthesis.jl")
include("../Utils/math_tools.jl")
include("leaf_energy_water_balance.jl")
include("photosynthesis_rates.jl")
include("rate_constants.jl")
include("respiration_models.jl")
include("fluorescence_models.jl")

end # module