# FT and NaN test
@testset "Photosynthesis --- FT consistency and not NaN" begin
    for FT in [Float32, Float64]
        envir  = PM.AirLayer{FT}();
        leaf   = PM.Leaf{FT}();
        mod_b  = PM.C3Bernacchi(FT);
        mod_3  = PM.C3CLM(FT);
        mod_4  = PM.C3CLM(FT);
        rand_T = rand(FT) + 298;

        leaf_b = deepcopy(leaf);
        leaf_3 = deepcopy(leaf);
        leaf_4 = deepcopy(leaf);

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

        for result in [ PM.arrhenius_correction(mod_3.KcT, rand_T),
                        PM.arrhenius_correction(mod_3.VcT, rand_T),
                        PM.photo_TD_from_set(mod_3.KcT, rand_T ),
                        PM.photo_TD_from_val(mod_3.VcT, leaf.Vcmax25 , rand_T)]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        PM.leaf_photo_from_pi!(mod_b, leaf_b, envir);
        PM.leaf_photo_from_pi!(mod_3, leaf_3, envir);
        PM.leaf_photo_from_pi!(mod_4, leaf_4, envir);
        recursive_FT_test(leaf_b, FT);
        recursive_FT_test(leaf_3, FT);
        recursive_FT_test(leaf_4, FT);
        recursive_NaN_test(leaf_b);
        recursive_NaN_test(leaf_3);
        recursive_NaN_test(leaf_4);

        PM.leaf_photo_from_glc!(mod_b, leaf_b, envir);
        PM.leaf_photo_from_glc!(mod_3, leaf_3, envir);
        PM.leaf_photo_from_glc!(mod_4, leaf_4, envir);
        recursive_FT_test(leaf_b, FT);
        recursive_FT_test(leaf_3, FT);
        recursive_FT_test(leaf_4, FT);
        recursive_NaN_test(leaf_b);
        recursive_NaN_test(leaf_3);
        recursive_NaN_test(leaf_4);

        for esm in [PM.ESMBallBerry{FT}(),
                    PM.ESMGentine{FT}(),
                    PM.ESMLeuning{FT}(),
                    PM.ESMMedlyn{FT}()]
            mod_b.Sto = esm;
            mod_3.Sto = esm;
            mod_4.Sto = esm;
            for res in [PM.empirical_gsw_from_model(esm, leaf_b, envir, FT(1)),
                        PM.empirical_gsw_from_model(esm, leaf_3, envir, FT(1)),
                        PM.empirical_gsw_from_model(esm, leaf_4, envir, FT(1))]
                recursive_FT_test(res, FT);
                recursive_NaN_test(res);
            end
        end

        for sm in [PM.ESMBallBerry{FT}(),
                   PM.ESMGentine{FT}(),
                   PM.ESMLeuning{FT}(),
                   PM.ESMMedlyn{FT}(),
                   PM.OSMEller(),
                   PM.OSMSperry(),
                   PM.OSMWang(),
                   PM.OSMWAP{FT}(),
                   PM.OSMWAPMod{FT}()]
            mod_b.Sto = sm;
            mod_3.Sto = sm;
            mod_4.Sto = sm;
            for res in [PM.envir_diff!(FT(20), mod_b, leaf_b, envir, mod_b.Sto),
                        PM.envir_diff!(FT(20), mod_3, leaf_3, envir, mod_3.Sto),
                        PM.envir_diff!(FT(20), mod_4, leaf_4, envir, mod_4.Sto)]
                recursive_FT_test(res, FT);
                recursive_NaN_test(res);
            end
            PM.leaf_photo_from_envir!(mod_b, leaf_b, envir, mod_b.Sto);
            PM.leaf_photo_from_envir!(mod_3, leaf_3, envir, mod_3.Sto);
            PM.leaf_photo_from_envir!(mod_4, leaf_4, envir, mod_4.Sto);
            PM.leaf_heat_flux!(leaf_b, envir);
            PM.leaf_heat_flux!(leaf_3, envir);
            PM.leaf_heat_flux!(leaf_4, envir);
            recursive_FT_test(leaf_b, FT);
            recursive_FT_test(leaf_3, FT);
            recursive_FT_test(leaf_4, FT);
            recursive_NaN_test(leaf_b);
            recursive_NaN_test(leaf_3);
            recursive_NaN_test(leaf_4);
        end
    end
end
