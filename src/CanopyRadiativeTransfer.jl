module CanopyRadiativeTransfer

using ClimaCache: BroadbandRadiation, BroadbandSLCanopy, BroadbandSoilAlbedo, HyperspectralMLCanopy, HyperspectralRadiation, HyperspectralSoilAlbedo, Leaves1D, Leaves2D, MonoMLGrassSPAC,
      MonoMLPalmSPAC, MonoMLTreeSPAC, Soil, SunSensorGeometry, VerhoefLIDF
using ConstrainedRootSolvers: ReduceStepMethodND, SolutionToleranceND, find_peak
using DocStringExtensions: METHODLIST
using LeafOptics: energy!, photon, photon!
using LinearAlgebra: mul!, pinv
using PkgUtility: K_STEFAN
using QuadGK: quadgk
using Statistics: mean
using UnPack: @unpack


const MODIS_BAND_1 = [ 620,  670];  # RED
const MODIS_BAND_2 = [ 841,  876];  # NIR
const MODIS_BAND_3 = [ 459,  479];  # BLUE
const MODIS_BAND_7 = [2105, 2155];  # SWIR

const OCO2_SIF_759 = [758.17, 759.20];  # SIF 757
const OCO2_SIF_770 = [769.62, 770.28];  # SIF 770
const OCO3_SIF_759 = [758.26, 759.28];  # SIF 757
const OCO3_SIF_770 = [769.67, 770.34];  # SIF 770

const TROPOMI_SIF_683 = [680, 685];     # SIF 683
const TROPOMI_SIF_747 = [735, 758];     # SIF 747
const TROPOMI_SIF_751 = [743, 758];     # SIF 751

const SOIL_ALBEDOS = [0.36 0.61 0.25 0.50;
                      0.34 0.57 0.23 0.46;
                      0.32 0.53 0.21 0.42;
                      0.31 0.51 0.20 0.40;
                      0.30 0.49 0.19 0.38;
                      0.29 0.48 0.18 0.36;
                      0.28 0.45 0.17 0.34;
                      0.27 0.43 0.16 0.32;
                      0.26 0.41 0.15 0.30;
                      0.25 0.39 0.14 0.28;
                      0.24 0.37 0.13 0.26;
                      0.23 0.35 0.12 0.24;
                      0.22 0.33 0.11 0.22;
                      0.20 0.31 0.10 0.20;
                      0.18 0.29 0.09 0.18;
                      0.16 0.27 0.08 0.16;
                      0.14 0.25 0.07 0.14;
                      0.12 0.23 0.06 0.12;
                      0.10 0.21 0.05 0.10;
                      0.08 0.16 0.04 0.08];

include("clumping.jl"         )
include("coefficients.jl"     )
include("fluorescence.jl"     )
include("geometry.jl"         )
include("inclination_angle.jl")
include("radiation.jl"        )
include("remote_sensing.jl"   )
include("soil.jl"             )


end # module
