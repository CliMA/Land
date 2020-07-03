using Test
using WaterPhysics

WP = WaterPhysics




# Test the variable FT recursively
function recursive_FT_test(para, FT)
    # if the type is AbstractFloat
    if typeof(para) <: AbstractFloat
        try
            @test typeof(para) == FT
        catch e
            println("The not NaN test failed for ", para, " and ", FT)
        end
    # if the type is array
    elseif typeof(para) <: AbstractArray
        if eltype(para) <: AbstractFloat
            try
                @test eltype(para) == FT
            catch e
                println("The not NaN test failed for ", para, " and ", FT)
            end
        else
            for ele in para
                recursive_FT_test(ele, FT)
            end
        end
    else
        # try if the parameter is a struct
        try
            for fn in fieldnames( typeof(para) )
                recursive_FT_test( getfield(para, fn), FT )
            end
        catch e
            println(typeof(para), "is not supprted by recursive_FT_test.")
        end
    end
end




# Test the variable NaN recursively
function recursive_NaN_test(para)
    # if the type is Number
    if typeof(para) <: Number
        try
            @test !isnan(para)
        catch e
            println("The not NaN test failed for", para)
        end
    # if the type is array
    elseif typeof(para) <: AbstractArray
        for ele in para
            recursive_NaN_test(ele)
        end
    else
        # try if the parameter is a struct
        try
            for fn in fieldnames( typeof(para) )
                recursive_NaN_test( getfield(para, fn) )
            end
        catch e
            println(typeof(para), "is not supprted by recursive_NaN_test.")
        end
    end
end




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
