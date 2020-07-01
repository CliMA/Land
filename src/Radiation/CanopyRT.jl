module CanopyRT

using DocStringExtensions
using MAT
using Parameters
using Polynomials
using QuadGK
using StaticArrays
using Statistics

using ..LandParameters
using ..MathTools

const K_BOLTZMANN = LandParameters.K_BOLTZMANN
const file_Opti   = joinpath(@__DIR__, "../../data/Optipar2017_ProspectD.mat")
const file_Sun    = joinpath(@__DIR__, "../../data/sun.mat")

# TODO remove the nWl, nAxi and etc. from types {}

# export types
export AbstractCanopyOpti,
       Canopy4RT,
       CanopyRadiation,
       SoilOpti,
       SolarAngles

# export public functions
export big_leaf_partition!,
       compute_canopy_geometry!,
       compute_canopy_matrices!,
       computeSIF_Fluxes!,
       create_canopy_optical,
       create_incoming_radiation,
       create_leaf_bio,
       create_opti_par,
       create_wl_para_set,
       derive_canopy_fluxes!,
       fluspect!,
       initialize_rt_module,
       simulate_short_wave!




include("types.jl"     )
include("bigleaf.jl"   )
include("simulation.jl")




###############################################################################
#
# Functions in development
# Test and document the functions before merging them to this file
#
###############################################################################
include("CanopyRT_in_development.jl")




end
