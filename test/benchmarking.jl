using BenchmarkTools
using WaterPhysics

function benchmark_WaterPhysics(FT)
    # define the variables
    rand_r = (rand(FT) + 20) * FT(1e-6);
    rand_T = rand(FT) + 298;
    rand_α = rand(FT) * 50;
    rand_Ψ = rand(FT) - 3;

    gas_air = TraceGasAir();
    gas_CO₂ = TraceGasCO₂();

    # benchmarking the functions
    println("\nUsing ", FT);
    println("\nBenchmarking capillary_pressure functions...");
    @btime capillary_pressure($rand_r, $rand_T);
    @btime capillary_pressure($rand_r, $rand_T, $rand_α);

    println("\nBenchmarking diffusive_coefficient functions...");
    @btime diffusive_coefficient($rand_T, $gas_CO₂, $gas_air);
    @btime relative_diffusive_coefficient($rand_T);

    println("\nBenchmarking latent_heat_vapor functions...");
    @btime latent_heat_vapor($rand_T);

    println("\nBenchmarking surface_tension functions...");
    @btime surface_tension($rand_T);
    @btime relative_surface_tension($rand_T);

    println("\nBenchmarking saturation_vapor_pressure functions...");
    @btime pressure_correction($rand_T, $rand_Ψ);
    @btime saturation_vapor_pressure($rand_T);
    @btime saturation_vapor_pressure($rand_T, $rand_Ψ);
    @btime saturation_vapor_pressure_slope($rand_T);
    @btime saturation_vapor_pressure_slope($rand_T, $rand_Ψ);

    println("\nBenchmarking viscosity functions...");
    @btime viscosity($rand_T);
    @btime relative_viscosity($rand_T);

    return nothing
end

benchmark_WaterPhysics(Float32);
benchmark_WaterPhysics(Float64);
