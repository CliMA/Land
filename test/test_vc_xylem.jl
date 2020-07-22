# FT and NaN test
println("\nTesting and benchmarking the xylem VC functions...")
@testset "Hydraulics --- xylem VC" begin
    for FT in [Float32, Float64]
        _vc1 = WeibullSingle{FT}();
        _vc2 = WeibullDual{FT}();
        _sh1 = BrooksCorey{FT}();
        _sh2 = VanGenuchten{FT}();
        _p   = FT(-0.5);
        _v   = FT(0.9);

        # Test xylem_p_crit
        f_st = FT(1)
        for result in [ xylem_p_crit(_vc1, f_st),
                        xylem_p_crit(_vc2, f_st) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        if benchmarking
            @btime xylem_p_crit($_vc1, $f_st);
            @btime xylem_p_crit($_vc2, $f_st);
        end

        # test the xylem_k_ratio
        for result in [ xylem_k_ratio(_vc1, _p, _v),
                        xylem_k_ratio(_vc2, _p, _v) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        if benchmarking
            @btime xylem_k_ratio($_vc1, $_p, $_v);
            @btime xylem_k_ratio($_vc2, $_p, $_v);
        end
    end
end
