# FT and NaN test
@testset "Leaf --- FT consistency and not NaN" begin
    for FT in [Float32, Float64]
        for data_set in [ LF.BLFixed(FT),
                      LF.LeafParams{FT}(),
                      LF.MeteoParams{FT}() ] 
            recursive_FT_test(data_set, FT)
            recursive_NaN_test(data_set, FT)
        end
    end
end




# Function test
@testset "Leaf --- Photosynthesis model" begin
    for FT in [Float32, Float64]
        C3_fluo = PM.FluorescenceFlexas(FT)
        C4_fluo = PM.FluorescenceFlexas(FT)
        C3_para = PM.C3CLM(FT)
        C4_para = PM.C4CLM(FT)
        C3_leaf = LF.LeafParams{FT}()
        C4_leaf = LF.LeafParams{FT}()

        LF.update_leaf_TD!(C3_para, C3_leaf)
        LF.update_leaf_TD!(C4_para, C4_leaf)
        LF.electron_transport_rate!(C3_para, C3_leaf, FT(500.0))
        LF.rubisco_limited_rate!(C3_para, C3_leaf)
        LF.rubisco_limited_rate!(C4_para, C4_leaf)
        LF.light_limited_rate!(C3_para, C3_leaf, FT(500.0))
        LF.light_limited_rate!(C4_para, C4_leaf, FT(500.0))
        LF.product_limited_rate!(C3_para, C3_leaf)
        LF.product_limited_rate!(C4_para, C4_leaf)
        LF.leaf_fluorescence!(C3_fluo, C3_leaf)
        LF.leaf_fluorescence!(C4_fluo, C4_leaf)

        @test true
    end
end
