using BenchmarkTools
using CanopyRadiation

function benchmark_CanopyRadiation(FT)
    println("\nBenchmarking the big_leaf_partition");
    _lai  = FT(2.0);
    _zen  = FT(30.0);
    _rall = FT(1000.0);
    @btime big_leaf_partition($_lai, $_zen, $_rall);
end

benchmark_CanopyRadiation(Float32);
benchmark_CanopyRadiation(Float64);
