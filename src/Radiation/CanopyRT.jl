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

# export types
export AbstractCanopyOpti,
       Canopy4RT,
       CanopyRadiation,
       SoilOpti,
       SolarAngles

# export public functions
export compute_canopy_geometry!,
       compute_canopy_matrices!,
       computeSIF_Fluxes!,
       create_canopy_optical,
       create_incoming_radiation,
       create_leaf_bio,
       create_opti_par,
       create_wl_para_set,
       derive_canopy_fluxes!,
       fluspect!,
       simulate_short_wave!




include("types.jl"     )
include("simulation.jl")




###############################################################################
#
# Functions in development
# Test and document the functions before merging them to this file
#
###############################################################################
include("CanopyRT_in_development.jl")




end
