using BenchmarkTools
using PlantHydraulics
using Test




include("recursive_test.jl")




benchmarking = false
include("test_struct.jl")
include("test_vc_soil.jl")
include("test_vc_xylem.jl")
include("test_leaf.jl")
include("test_legacy.jl")
include("test_temperature.jl")
include("test_root.jl")
include("test_pressure.jl")
include("test_plant.jl")

benchmarking = true
