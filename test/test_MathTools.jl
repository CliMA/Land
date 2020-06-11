# FT and NaN test
@testset "MathTools --- FT consistency and not NaN" begin
    for FT in [Float32, Float64]
        for result in [ MT.dcum(rand(FT)/2, rand(FT)/2, rand(FT)*90),
                        MT.dcum(rand(FT)/2, rand(FT)/2, rand(FT)*90),
                        MT.dladgen(rand(FT)/2, rand(FT)/2, rand(FT,10,2)),
                        MT.e2phot(rand(FT, 100), rand(FT,100)),
                        MT.fast∫(rand(FT, 100), rand(FT,100)),
                        MT.psofunction(rand(FT), rand(FT), rand(FT), rand(FT), rand(FT), rand(FT), rand(FT)),
                        MT.lower_quadratic(rand(FT), rand(FT)+2, rand(FT)),
                        MT.quadratic(rand(FT), rand(FT)+2, rand(FT)),
                        MT.volscatt(rand(FT), rand(FT), rand(FT), rand(FT)) ]
            recursive_FT_test(result, FT)
            recursive_NaN_test(result, FT)
        end
    end
end
