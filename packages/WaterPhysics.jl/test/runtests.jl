using PkgUtility
using Test
using WaterPhysics


@testset verbose = true "WaterPhysics" begin
    @testset "Capillary Pressure" begin
        for FT in [Float32, Float64]
            rand_r = (rand(FT) + 20) * FT(1e-6);
            rand_T = rand(FT) + 298;
            rand_α = rand(FT) * 50;
            for result in [ WaterPhysics.capillary_pressure(rand_r, rand_T),
                            WaterPhysics.capillary_pressure(rand_r, rand_T, rand_α) ]
                @test PkgUtility.FT_test(result, FT);
                @test PkgUtility.NaN_test(result);
            end;
        end;
    end;

    @testset "Diffusive Coefficient" begin
        for FT in [Float32, Float64]
            rand_T  = rand(FT) + 298;
            gas_air = WaterPhysics.TraceGasAir{FT}();
            gas_CO₂ = WaterPhysics.TraceGasCO₂{FT}();
            liq_H₂O = WaterPhysics.TraceLiquidH₂O{FT}();
            for result in [ WaterPhysics.diffusive_coefficient(rand_T, gas_CO₂, gas_air),
                            WaterPhysics.diffusive_coefficient(rand_T, gas_CO₂, liq_H₂O),
                            WaterPhysics.relative_diffusive_coefficient(rand_T),
                            WaterPhysics.relative_diffusive_coefficient(rand_T, gas_CO₂, liq_H₂O) ]
                @test PkgUtility.FT_test(result, FT);
                @test PkgUtility.NaN_test(result);
            end;
            @test WaterPhysics.relative_diffusive_coefficient(FT(298.15)) ≈ 1;
        end;
    end;

    @testset "Latent Heat of Vapor" begin
        for FT in [Float32, Float64]
            rand_T = rand(FT) + 298;
            for result in [ WaterPhysics.latent_heat_vapor(rand_T) ]
                @test PkgUtility.FT_test(result, FT);
                @test PkgUtility.NaN_test(result);
            end
        end;
    end;

    @testset "Surface Tension" begin
        for FT in [Float32, Float64]
            rand_T = rand(FT) + 298;
            for result in [ WaterPhysics.surface_tension(rand_T),
                            WaterPhysics.relative_surface_tension(rand_T) ]
                @test PkgUtility.FT_test(result, FT);
                @test PkgUtility.NaN_test(result);
            end;
            @test WaterPhysics.relative_surface_tension(FT(298.15)) ≈ 1;
        end;
    end;

    @testset "Saturated Vapor Pressure" begin
        for FT in [Float32, Float64]
            rand_T = rand(FT) + 298;
            rand_Ψ = rand(FT) - 3;
            for result in [ WaterPhysics.pressure_correction(rand_T, rand_Ψ),
                            WaterPhysics.saturation_vapor_pressure(rand_T),
                            WaterPhysics.saturation_vapor_pressure(rand_T, rand_Ψ),
                            WaterPhysics.saturation_vapor_pressure_slope(rand_T),
                            WaterPhysics.saturation_vapor_pressure_slope(rand_T, rand_Ψ) ]
                @test PkgUtility.FT_test(result, FT);
                @test PkgUtility.NaN_test(result);
            end;
        end;
    end;

    @testset "Viscosity" begin
        for FT in [Float32, Float64]
            rand_T = rand(FT) + 298;
            for result in [ WaterPhysics.viscosity(rand_T),
                            WaterPhysics.relative_viscosity(rand_T) ]
                @test PkgUtility.FT_test(result, FT);
                @test PkgUtility.NaN_test(result);
            end;
            @test WaterPhysics.relative_viscosity(FT(298.15)) ≈ 1;
        end;
    end;
end;
