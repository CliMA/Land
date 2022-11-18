@testset verbose = true "EarthSurface Test" begin
    for FT in [Float32, Float64]
        @test 0 <= EarthSurface.solar_zenith_angle(FT(20), FT(12.5)) <= 90;
        @test 0 <= EarthSurface.solar_zenith_angle(FT(20), FT(12), FT(12), FT(45)) <= 90;
    end;
end;
