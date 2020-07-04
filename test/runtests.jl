using Photosynthesis
using Test

PM = Photosynthesis




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
@testset "Photosynthesis --- FT consistency and not NaN" begin
    for FT in [Float32, Float64]
        envir  = PM.AirLayer{FT}();
        leaf   = PM.Leaf{FT}();
        mod_b  = PM.C3Bernacchi(FT);
        mod_3  = PM.C3CLM(FT);
        mod_4  = PM.C4CLM(FT);
        rand_T = rand(FT) + 298;

        leaf_b = deepcopy(leaf);
        leaf_3 = deepcopy(leaf);
        leaf_4 = deepcopy(leaf);

        # Test types
        for data_set in [ envir,
                          leaf,
                          mod_b,
                          mod_3,
                          mod_4,
                          PM.JmaxTDLeuning(FT),
                          PM.KpepTDBoyd(FT),
                          PM.VcmaxTDLeuning(FT),
                          PM.VtoRCollatz(FT)]
            recursive_FT_test(data_set, FT);
            recursive_NaN_test(data_set);
        end

        # Test temperature dependence
        for result in [ PM.arrhenius_correction(mod_3.KcT, rand_T),
                        PM.arrhenius_correction(mod_3.VcT, rand_T),
                        PM.photo_TD_from_set(mod_3.KcT, rand_T ),
                        PM.photo_TD_from_val(mod_3.VcT, leaf.Vcmax25 , rand_T)]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        # test photosynthesis model using p_i
        PM.leaf_photo_from_pi!(mod_b, leaf_b, envir);
        PM.leaf_photo_from_pi!(mod_3, leaf_3, envir);
        PM.leaf_photo_from_pi!(mod_4, leaf_4, envir);
        recursive_FT_test(leaf_b, FT);
        recursive_FT_test(leaf_3, FT);
        recursive_FT_test(leaf_4, FT);
        recursive_NaN_test(leaf_b);
        recursive_NaN_test(leaf_3);
        recursive_NaN_test(leaf_4);

        # test photosynthesis model using glc
        PM.leaf_photo_from_glc!(mod_b, leaf_b, envir);
        PM.leaf_photo_from_glc!(mod_3, leaf_3, envir);
        PM.leaf_photo_from_glc!(mod_4, leaf_4, envir);
        recursive_FT_test(leaf_b, FT);
        recursive_FT_test(leaf_3, FT);
        recursive_FT_test(leaf_4, FT);
        recursive_NaN_test(leaf_b);
        recursive_NaN_test(leaf_3);
        recursive_NaN_test(leaf_4);

        # test photosynthesis model for 1D Array of glc
        glcs = rand(FT, 100) ./ 20 .+ FT(0.1);
        pars = rand(FT, 100) .+ 1000;
        Jps  = PM.leaf_ETR_pot_APAR(leaf, pars);
        Js3  = PM.leaf_ETR_Jps(mod_3, leaf, Jps);
        Js4  = PM.leaf_ETR_Jps(mod_4, leaf, Jps);
        Acs3 = PM.rubisco_limited_an_glc(mod_3, leaf, envir, glcs);
        Ajs3 = PM.light_limited_an_glc(mod_3, leaf, envir, glcs, Js3);
        Aps3 = PM.product_limited_an_glc(mod_3, leaf, envir, glcs);
        Acs4 = PM.rubisco_limited_an_glc(mod_4, leaf, envir, glcs);
        Ajs4 = PM.light_limited_an_glc(mod_4, leaf, envir, glcs, Js4);
        Aps4 = PM.product_limited_an_glc(mod_4, leaf, envir, glcs);
        for result in [ Jps, Js3, Js4, Acs3, Ajs3, Aps3, Acs4, Ajs4, Aps4]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end
    end
end
