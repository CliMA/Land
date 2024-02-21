@testset "CLM5 Mode" begin
    FT = Float64;
    wls = WaveLengths{FT}(maxwlPAR = 750);
    node = SoilPlantAirContinuum.SPACMono{FT}(wl_set = wls);
    for leaf in node.leaves_rt
        leaf.prescribe = true;
    end;
    node.soil_opt.hyperspectral = false;
    SoilPlantAirContinuum.initialize_spac_canopy!(node);
    @test true;
end;
