module WaterPhysics

using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using UnPack: @unpack

using EmeraldConstants: CP_L, CP_V, GAS_R, K_BOLTZMANN, LH_V₀, PRESS_TRIPLE, R_V, T₂₅, T_TRIPLE, V_H₂O


include("water.jl")


end # module
