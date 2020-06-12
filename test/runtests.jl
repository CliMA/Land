using Statistics
using Test

ENV["JULIA_LOG_LEVEL"] = "WARN"




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
function recursive_NaN_test(para, FT)
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
            recursive_NaN_test(ele, FT)
        end
    else
        # try if the parameter is a struct
        try
            for fn in fieldnames( typeof(para) )
                recursive_NaN_test( getfield(para, fn), FT )
            end
        catch e
            println(typeof(para), "is not supprted by recursive_NaN_test.")
        end
    end
end




# When moving any function from developing.jl to stable files,
# add corresponding test and docs
using Land

MT   = Land.MathTools
WP   = Land.WaterPhysics
PM   = Land.Photosynthesis
LF   = Land.Leaf
PT   = Land.Plant
RT   = Land.CanopyRT
SPAC = Land.SPAC

include("test_MathTools.jl"     )
include("test_WaterPhysics.jl"  )
include("test_Photosynthesis.jl")
include("test_Leaf.jl"          )
include("test_Plant.jl"         )
include("test_CanopyRT.jl"      )
include("test_SPAC.jl"          )

println("All tests finished.")
