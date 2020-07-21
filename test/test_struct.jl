# FT and NaN test for the structs
println("Testing FT and NaN for the structures...");
@testset "Photosynthesis --- struct" begin
    for FT in [Float32, Float64]
        for data_set in [ KcTDBernacchi(FT),
                          VpmaxTDBoyd(FT),
                          C3CLM(FT),
                          C4CLM(FT),
                          AirLayer{FT}(),
                          Leaf{FT}() ]
            recursive_FT_test(data_set, FT);
            recursive_NaN_test(data_set);
        end
    end
end
