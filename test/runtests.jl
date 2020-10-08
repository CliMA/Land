using Test
using WaterPhysics




include("recursive_test.jl");




@testset "WaterPhysics --- functions" begin
    for FT in [Float32, Float64]
        rand_r = (rand(FT) + 20) * FT(1e-6);
        rand_T = rand(FT) + 298;
        rand_α = rand(FT) * 50;
        rand_Ψ = rand(FT) - 3;

        gas_air = TraceGasAir();
        gas_CO₂ = TraceGasCO₂();

        # capillary_pressure
        for result in [ capillary_pressure(rand_r, rand_T),
                        capillary_pressure(rand_r, rand_T, rand_α) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        # diffusive_coefficient
        for result in [ diffusive_coefficient(rand_T, gas_CO₂, gas_air),
                        relative_diffusive_coefficient(rand_T) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end
        @test relative_diffusive_coefficient(FT(298.15)) ≈ 1;

        # latent_heat_vapor
        for result in [ latent_heat_vapor(rand_T) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        # surface_tension and relative_surface_tension
        for result in [ surface_tension(rand_T),
                        relative_surface_tension(rand_T) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end
        @test relative_surface_tension(FT(298.15)) ≈ 1;

        # saturation_vapor_pressure
        for result in [ pressure_correction(rand_T, rand_Ψ),
                        saturation_vapor_pressure(rand_T),
                        saturation_vapor_pressure(rand_T, rand_Ψ),
                        saturation_vapor_pressure_slope(rand_T),
                        saturation_vapor_pressure_slope(rand_T, rand_Ψ) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        # viscosity and relative_viscosity
        for result in [ viscosity(rand_T),
                        relative_viscosity(rand_T) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end
        @test relative_viscosity(FT(298.15)) ≈ 1;
    end
end
