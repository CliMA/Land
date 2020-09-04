using BenchmarkTools
using Photosynthesis

function benchmark_Photosynthesis(FT)
    # define variables
    c3_set   = C3CLM(FT);
    c4_set   = C4CLM(FT);
    leaf_3   = Leaf{FT}();
    leaf_4   = Leaf{FT}();
    envir    = AirLayer{FT}();
    fluo_set = c3_set.Flu;
    T        = rand(FT) + 298;
    glc      = FT(0.1);
    p_i      = rand(FT) + 20;

    # benchmarking the functions
    println("\nUsing ", FT);
    println("\nBenchmarking leaf_temperature_dependence! functions...");
    @btime leaf_temperature_dependence!($c3_set, $leaf_3, $envir, $T);
    @btime leaf_temperature_dependence!($c4_set, $leaf_4, $envir, $T);

    println("\nBenchmarking rubisco_limited_rate*! functions...");
    @btime rubisco_limited_rate!($c3_set, $leaf_3);
    @btime rubisco_limited_rate!($c4_set, $leaf_4);
    @btime rubisco_limited_rate_glc!($c3_set, $leaf_3, $envir);

    println("\nBenchmarking leaf_ETR! functions...");
    @btime leaf_ETR!($c3_set, $leaf_3);
    @btime leaf_ETR!($c4_set, $leaf_4);

    println("\nBenchmarking light_limited_rate! functions...");
    @btime light_limited_rate!($c3_set, $leaf_3);
    @btime light_limited_rate!($c4_set, $leaf_4);
    @btime light_limited_rate_glc!($c3_set, $leaf_3, $envir);

    println("\nBenchmarking product_limited_rate! functions...");
    @btime product_limited_rate!($c3_set, $leaf_3);
    @btime product_limited_rate!($c4_set, $leaf_4);
    @btime product_limited_rate_glc!($c4_set, $leaf_4, $envir);

    println("\nBenchmarking leaf_fluorescence! functions...");
    @btime leaf_fluorescence!($fluo_set, $leaf_3);

    println("\nBenchmarking leaf_photo_from_glc! functions...");
    @btime leaf_photo_from_glc!($c3_set, $leaf_3, $envir, $glc);
    @btime leaf_photo_from_glc!($c4_set, $leaf_4, $envir, $glc);

    println("\nBenchmarking leaf_photo_from_pi! functions...");
    @btime leaf_photo_from_pi!($c3_set, $leaf_3, $p_i);
    @btime leaf_photo_from_pi!($c4_set, $leaf_4, $p_i);
end

benchmark_Photosynthesis(Float32);
benchmark_Photosynthesis(Float64);
