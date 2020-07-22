# FT and NaN test
println("\nTesting and benchmarking the temperature functions...")
@testset "Hydraulics --- temperature functions" begin
    for FT in [Float32, Float64]
        leaf  = LeafHydraulics{FT}();
        root  = RootHydraulics{FT}();
        treet = TreeSimple{FT}();
        T     = rand(FT) + 298;

        # test the temperature functions
        vc_temperature_effects!(treet);
        vc_temperature_effects!(leaf, T);
        vc_temperature_effects!(root, T);

        for dataset in [leaf, root, treet]
            recursive_FT_test(dataset, FT);
            recursive_NaN_test(dataset);
        end

        if benchmarking
            @btime vc_temperature_effects!($treet);
            @btime vc_temperature_effects!($leaf, $T);
            @btime vc_temperature_effects!($root, $T);
        end
    end
end
