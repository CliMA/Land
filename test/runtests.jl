using BenchmarkTools
using Photosynthesis
using Test




include("recursive_test.jl")




# FT and NaN test for the structs
@testset "Photosynthesis --- structures" begin
    for FT in [Float32, Float64]
        for data_set in [ KcTDBernacchi(FT),
                          VpmaxTDBoyd(FT),
                          C3CLM(FT),
                          C4CLM(FT),
                          AirLayer{FT}(),
                          Leaf{FT}(),
                          KoTDBernacchi(FT),
                          RespirationTDBernacchi(FT),
                          VcmaxTDBernacchi(FT),
                          VomaxTDBernacchi(FT),
                          ΓStarTDBernacchi(FT),
                          KpepTDBoyd(FT),
                          JmaxTDLeuning(FT),
                          VcmaxTDLeuning(FT),
                          JmaxTDBernacchi(FT),
                          VtoRCollatz(FT),
                          C3Bernacchi(FT) ]
            recursive_FT_test(data_set, FT);
            recursive_NaN_test(data_set);
        end
    end
end




# FT and NaN test for the structs
@testset "Photosynthesis --- functions" begin
    for FT in [Float32, Float64]
        c3_set   = C3CLM(FT);
        c4_set   = C4CLM(FT);
        leaf_3   = Leaf{FT}();
        leaf_4   = Leaf{FT}();
        td_q10   = Q10TD{FT}(1, 273.15, 1.7);
        envir    = AirLayer{FT}();
        fluo_set = c3_set.Flu;
        T        = rand(FT) + 298;
        glc      = FT(0.1);
        p_i      = rand(FT) + 20;

        # temperature corrections
        photo_TD_from_set(td_q10, T);
        leaf_temperature_dependence!(c3_set, leaf_3, envir, T);
        leaf_temperature_dependence!(c4_set, leaf_4, envir, T);

        # rubisco limited rates
        rubisco_limited_rate!(c3_set, leaf_3);
        rubisco_limited_rate!(c4_set, leaf_4);
        recursive_NaN_test(leaf_3);
        recursive_NaN_test(leaf_4);
        rubisco_limited_rate_glc!(c3_set, leaf_3, envir);
        recursive_NaN_test(leaf_3);

        # light limited rates
        leaf_ETR!(c3_set, leaf_3);
        leaf_ETR!(c4_set, leaf_4);
        light_limited_rate!(c3_set, leaf_3);
        light_limited_rate!(c4_set, leaf_4);
        recursive_NaN_test(leaf_3);
        recursive_NaN_test(leaf_4);
        light_limited_rate_glc!(c3_set, leaf_3, envir);
        recursive_NaN_test(leaf_3);

        # product limited rates
        product_limited_rate!(c3_set, leaf_3);
        product_limited_rate!(c4_set, leaf_4);
        recursive_NaN_test(leaf_3);
        recursive_NaN_test(leaf_4);
        product_limited_rate_glc!(c4_set, leaf_4, envir);
        recursive_NaN_test(leaf_4);

        # fluorescence
        leaf_photo_from_glc!(c3_set, leaf_3, envir);
        leaf_fluorescence!(fluo_set, leaf_3);
        recursive_NaN_test(leaf_3);

        # leaf photo from glc
        leaf_photo_from_glc!(c3_set, leaf_3, envir, glc);
        leaf_photo_from_glc!(c4_set, leaf_4, envir, glc);
        recursive_NaN_test(leaf_3);
        recursive_NaN_test(leaf_4);

        # leaf photo from p_i
        leaf_photo_from_pi!(c3_set, leaf_3, p_i);
        leaf_photo_from_pi!(c4_set, leaf_4, p_i);
        recursive_NaN_test(leaf_3);
        recursive_NaN_test(leaf_4);
    end
end
