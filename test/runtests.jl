using Photosynthesis
using PlantHydraulics
using StomataModels
using Test

PM = Photosynthesis
SM = StomataModels




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




# FT and NaN test
@testset "StomtaModels --- FT consistency and not NaN" begin
    for FT in [Float32, Float64]
        envir  = PM.AirLayer{FT}();
        leaf_3 = PM.Leaf{FT}();
        leaf_4 = PM.Leaf{FT}();
        mod_3  = PM.C3CLM(FT);
        mod_4  = PM.C4CLM(FT);
        rand_T = rand(FT) + 298;
        lv_3   = SM.Leaves{FT}(n_leaf=2);
        lv_4   = SM.Leaves{FT}(n_leaf=2);
        esm_1  = SM.ESMBallBerry{FT}();
        esm_2  = SM.ESMGentine{FT}();
        esm_3  = SM.ESMLeuning{FT}();
        esm_4  = SM.ESMMedlyn{FT}();
        osm_1  = SM.OSMEller{FT}();
        osm_2  = SM.OSMSperry{FT}();
        osm_3  = SM.OSMWang{FT}();
        osm_4  = SM.OSMWAP{FT}();
        osm_5  = SM.OSMWAPMod{FT}();
        hs     = LeafHydraulics{FT}();

        # test the structures
        for result in [ lv_3, esm_1, esm_2, esm_3, esm_4, osm_4, osm_5 ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        # test the refresh functions
        SM.update_leaf_TP!(mod_3, lv_3, hs, envir);
        SM.update_leaf_TP!(mod_4, lv_4, hs, envir);
        SM.update_leaf_AK!(mod_3, lv_3, hs, envir);
        SM.update_leaf_AK!(mod_4, lv_4, hs, envir);
        SM.update_leaf_from_glc!(mod_3, lv_3, envir, 1, FT(0.1));
        SM.update_leaf_from_glc!(mod_4, lv_4, envir, 1, FT(0.1));
        SM.update_leaf_from_gsw!(mod_3, lv_3, envir, 1, FT(0.05));
        SM.update_leaf_from_gsw!(mod_4, lv_4, envir, 1, FT(0.05));
        lv_3.g_sw[2] = 0;
        lv_4.g_sw[2] = 0;
        SM.leaf_gsw_control!(mod_3, lv_3, envir, 2);
        SM.leaf_gsw_control!(mod_4, lv_4, envir, 2);
        recursive_NaN_test(lv_3);
        recursive_NaN_test(lv_4);

        # test the empirical model formulations
        for result in [ SM.empirical_gsw_from_model(esm_1, leaf_3, envir, FT(1)),
                        SM.empirical_gsw_from_model(esm_2, leaf_3, envir, FT(1)),
                        SM.empirical_gsw_from_model(esm_3, leaf_3, envir, FT(1)),
                        SM.empirical_gsw_from_model(esm_4, leaf_3, envir, FT(1)),
                        SM.empirical_gsw_from_model(esm_1, lv_3, envir, FT(1)),
                        SM.empirical_gsw_from_model(esm_2, lv_3, envir, FT(1)),
                        SM.empirical_gsw_from_model(esm_3, lv_3, envir, FT(1)),
                        SM.empirical_gsw_from_model(esm_4, lv_3, envir, FT(1)),
                        SM.empirical_gsw_from_model(esm_1, lv_4, envir, FT(1)),
                        SM.empirical_gsw_from_model(esm_2, lv_4, envir, FT(1)),
                        SM.empirical_gsw_from_model(esm_3, lv_4, envir, FT(1)),
                        SM.empirical_gsw_from_model(esm_4, lv_4, envir, FT(1)),
                        SM.empirical_gsw_from_model(esm_1, lv_3, envir, FT(1), 1),
                        SM.empirical_gsw_from_model(esm_2, lv_3, envir, FT(1), 1),
                        SM.empirical_gsw_from_model(esm_3, lv_3, envir, FT(1), 1),
                        SM.empirical_gsw_from_model(esm_4, lv_3, envir, FT(1), 1),
                        SM.empirical_gsw_from_model(esm_1, lv_4, envir, FT(1), 1),
                        SM.empirical_gsw_from_model(esm_2, lv_4, envir, FT(1), 1),
                        SM.empirical_gsw_from_model(esm_3, lv_4, envir, FT(1), 1),
                        SM.empirical_gsw_from_model(esm_4, lv_4, envir, FT(1), 1) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        # test the solution functions
        for (mod,lv) in zip([mod_3, mod_3], [lv_3, lv_4])
            for sm in [esm_1, esm_2, esm_3, esm_4, osm_1, osm_2, osm_3, osm_4, osm_5]
                result = SM.envir_diff!(FT(0.1), mod, lv, hs, envir, sm, 1);
                recursive_FT_test(result, FT);
                recursive_NaN_test(result);
            end
        end

        # test the stomata solutions
        for (mod,lv) in zip([mod_3, mod_3], [lv_3, lv_4])
            for sm in [esm_1, esm_2, esm_3, esm_4, osm_1, osm_2, osm_3, osm_4, osm_5]
                SM.leaf_photo_from_envir!(mod, lv, hs, envir, sm, 1);
                recursive_NaN_test(lv);
            end
        end
        for (mod,lv) in zip([mod_3, mod_3], [lv_3, lv_4])
            for sm in [esm_1, esm_2, esm_3, esm_4, osm_1, osm_2, osm_3, osm_4, osm_5]
                SM.leaf_photo_from_envir!(mod, lv, hs, envir, sm);
                recursive_NaN_test(lv);
            end
        end
    end
end
