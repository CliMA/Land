# FT and NaN test
@testset "WaterPhysics --- FT consistency and not NaN" begin
    for FT in [Float32, Float64]
        rand_T  = rand(FT) + 298
        rand_Tl = rand(FT,10) .+ 298
        for result in [ WP.latent_heat_vapor(rand_T),
                        WP.saturation_vapor_pressure(rand_T),
                        WP.saturation_vapor_pressure_slope(rand_T),
                        WP.relative_surface_tension(rand_T),
                        WP.surface_tension(rand_T),
                        WP.relative_viscosity(rand_T),
                        WP.viscosity(rand_T),
                        WP.relative_diffusive_coefficient(rand_T),
                        WP.latent_heat_vapor(rand_Tl),
                        WP.saturation_vapor_pressure(rand_Tl),
                        WP.saturation_vapor_pressure_slope(rand_Tl),
                        WP.relative_surface_tension(rand_Tl),
                        WP.surface_tension(rand_Tl),
                        WP.relative_viscosity(rand_Tl),
                        WP.viscosity(rand_Tl),
                        WP.relative_diffusive_coefficient(rand_Tl) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end
    end
end
