# FT and NaN test
println("\nTesting and benchmarking the root-related functions...")
@testset "Hydraulics --- root-related functions" begin
    for FT in [Float32, Float64]
        root  = RootHydraulics{FT}();
        grass = create_grass_like_hs(FT(-2.1), FT(0.5), collect(FT,0:-0.5:-3.0), collect(FT,0:1:20));

        # test the root q functions
        _p1 = FT(0)
        _p2 = FT(-0.5)
        _p3 = FT(-1.0)
        for result in [ root_q_from_pressure(root, _p1),
                        root_q_from_pressure(root, _p2),
                        root_q_from_pressure(root, _p3) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        if benchmarking
            @btime root_q_from_pressure($root, $_p2);
        end

        # test the root qs from q functions
        _f0 = FT(0);
        _f1 = FT(0.1);
        _f2 = FT(0.5);
        _ks = zeros(FT, 5);
        _ps = zeros(FT, 5);
        _qs = zeros(FT, 5);
        for _f in [_f0, _f1, _f2]
            roots_flow!(grass.roots, _ks, _ps, _qs, _f);
            recursive_FT_test(grass, FT);
            recursive_NaN_test(grass);
        end

        # test the recalculate_roots_flow function
        _ks = zeros(FT, 5);
        _ps = zeros(FT, 5);
        _qs = zeros(FT, 5);
        recalculate_roots_flow!(grass.roots, _ks, _ps, _qs, FT(0.5));
        recursive_NaN_test(grass);

        if benchmarking
            grass = create_grass_like_hs(FT(-2.1), FT(0.5), collect(FT,0:-0.5:-3.0), collect(FT,0:1:20));
            @btime recalculate_roots_flow!($(grass.roots), $_ks, $_ps, $_qs, $_f2);
            @btime roots_flow!($(grass.roots), $_ks, $_ps, $_qs, $_f2);
        end
    end
end
