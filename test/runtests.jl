using ClimaCache
using SoilPlantAirContinuum
using Test


@testset verbose = true "SoilPlantAirContinuum Test" begin
    @testset "Core Function" begin
        for FT in [Float32, Float64]
            for spac in [ClimaCache.MonoMLGrassSPAC{FT}(),
                         ClimaCache.MonoMLPalmSPAC{FT}(),
                         ClimaCache.MonoMLTreeSPAC{FT}()]
                SoilPlantAirContinuum.initialize!(spac);
                SoilPlantAirContinuum.soil_plant_air_continuum!(spac, FT(1));
                @test true;
            end;
        end;
    end;
end;


#=

using ClimaCache, PlantHydraulics, SoilPlantAirContinuum;
FT = Float64;
spac = ClimaCache.MonoMLTreeSPAC{FT}();

for _leaf in spac.LEAVES
    _leaf.PSM.v_cmax25 = 50;
    _leaf.PSM.r_d25 = 0.75;
    _leaf.HS.AREA = 75;
end

SoilPlantAirContinuum.initialize!(spac);

SoilPlantAirContinuum.soil_plant_air_continuum!(spac, FT(10));

[soil.θ for soil in spac.SOIL.LAYERS]' * abs.(diff(spac.SOIL.ZS))
mean(spac.LEAVES[1].a_net_sunlit)
[leaf.t for leaf in spac.LEAVES]

var = [_leaf.HS.p_leaf for _leaf in spac.LEAVES]
@info var;

spac.ROOTS[2].HS.p_ups

=#
