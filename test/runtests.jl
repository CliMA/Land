using BenchmarkTools
using Photosynthesis
using PlantHydraulics
using StomataModels
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
include("test_struct.jl"   )
include("test_empirical.jl")

benchmarking = true
include("test_solution.jl")

#=
# FT and NaN test
@testset "StomtaModels --- FT consistency and not NaN" begin
    for FT in [Float32, Float64]
        envir  = AirLayer{FT}();
        leaf_3 = Leaf{FT}();
        leaf_4 = Leaf{FT}();
        mod_3  = C3CLM(FT);
        mod_4  = C4CLM(FT);
        rand_T = rand(FT) + 298;
        lv_3   = CanopyLayer{FT}(n_leaf=2);
        lv_4   = CanopyLayer{FT}(n_leaf=2);
        esm_1  = ESMBallBerry{FT}();
        esm_2  = ESMGentine{FT}();
        esm_3  = ESMLeuning{FT}();
        esm_4  = ESMMedlyn{FT}();
        osm_1  = OSMEller{FT}();
        osm_2  = OSMSperry{FT}();
        osm_3  = OSMWang{FT}();
        osm_4  = OSMWAP{FT}();
        osm_5  = OSMWAPMod{FT}();
        hs     = LeafHydraulics{FT}();

        # test the refresh functions
        update_leaf_TP!(mod_3, lv_3, hs, envir);
        update_leaf_TP!(mod_4, lv_4, hs, envir);
        update_leaf_AK!(mod_3, lv_3, hs, envir);
        update_leaf_AK!(mod_4, lv_4, hs, envir);
        update_leaf_from_glc!(mod_3, lv_3, envir, 1, FT(0.1));
        update_leaf_from_glc!(mod_4, lv_4, envir, 1, FT(0.1));
        update_leaf_from_gsw!(mod_3, lv_3, envir, 1, FT(0.05));
        update_leaf_from_gsw!(mod_4, lv_4, envir, 1, FT(0.05));
        lv_3.g_sw[2] = 0;
        lv_4.g_sw[2] = 0;
        leaf_gsw_control!(mod_3, lv_3, envir, 2);
        leaf_gsw_control!(mod_4, lv_4, envir, 2);
        recursive_NaN_test(lv_3);
        recursive_NaN_test(lv_4);
    end
end
=#