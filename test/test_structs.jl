# FT and NaN tests
# TODO test the MArrays as well
@testset "CanopyRT --- FT and NaN test" begin
    for FT in [Float32, Float64]
        for data_set in [ Canopy4RT{FT}(nLayer=20, LAI=FT(3)),
                          create_canopy_opticals(FT, 10, 10, 10, 10),
                          CanopyRads{FT}(nWL=10, nWLf=10, nLayer=10, nAzi=10, nIncl=9),
                          create_incoming_radiation(FT[600,620,640]),
                          create_leaf_bios(FT, 10, 10, 10),
                          create_leaf_opticals(FT[600,620,640]),
                          SoilOpticals{FT}(FT[600,610,620], FT[0.1,0.1,0.1], FT[0.1,0.1,0.1], FT(280.0)),
                          SolarAngles{FT}(),
                          WaveLengths{FT}() ]
            recursive_FT_test(data_set, FT)
            recursive_NaN_test(data_set)
        end
    end
end
