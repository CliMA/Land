using BenchmarkTools
using PlantHydraulics
using Test




# Test the variable FT recursively
function recursive_FT_test(para, FT)
    # if the type is AbstractFloat
    if typeof(para) <: AbstractFloat
        try
            @test typeof(para) == FT
        catch e
            println("The not NaN test failed for ", para, " and ", FT)
        end
    # if the type is array
    elseif typeof(para) <: AbstractArray
        if eltype(para) <: AbstractFloat
            try
                @test eltype(para) == FT
            catch e
                println("The not NaN test failed for ", para, " and ", FT)
            end
        else
            for ele in para
                recursive_FT_test(ele, FT)
            end
        end
    else
        # try if the parameter is a struct
        try
            for fn in fieldnames( typeof(para) )
                recursive_FT_test( getfield(para, fn), FT )
            end
        catch e
            println(typeof(para), "is not supprted by recursive_FT_test.")
        end
    end
end




# Test the variable NaN recursively
function recursive_NaN_test(para)
    # if the type is Number
    if typeof(para) <: Number
        try
            @test !isnan(para)
        catch e
            println("The not NaN test failed for", para)
        end
    # if the type is array
    elseif typeof(para) <: AbstractArray
        for ele in para
            recursive_NaN_test(ele)
        end
    else
        # try if the parameter is a struct
        try
            for fn in fieldnames( typeof(para) )
                recursive_NaN_test( getfield(para, fn) )
            end
        catch e
            println(typeof(para), "is not supprted by recursive_NaN_test.")
        end
    end
end




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
