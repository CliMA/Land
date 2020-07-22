# FT and NaN test
println("\nTesting and benchmarking the legacy functions...")
@testset "Hydraulics --- legacy" begin
    for FT in [Float32, Float64]
        _ps   = FT(-0.1);
        _ft   = FT(1);
        _fsl  = FT(0.8);
        _fsh  = FT(0.2);
        _rsl  = FT(0.5);

        # Test the struct
        treet = TreeSimple{FT}();
        hydraulic_p_profile!(treet, _ps, _ft);
        recursive_FT_test(treet, FT);
        recursive_NaN_test(treet);

        treet = TreeSimple{FT}();
        hydraulic_p_profile!(treet, _ps, _fsl, _fsh, _rsl);
        recursive_FT_test(treet, FT);
        recursive_NaN_test(treet);

        inititialize_legacy!(treet);
        recursive_FT_test(treet, FT);
        recursive_NaN_test(treet);

        if benchmarking
            @btime hydraulic_p_profile!($treet, $_ps, $_ft);
            @btime hydraulic_p_profile!($treet, $_ps, $_fsl, $_fsh, $_rsl);
            @btime inititialize_legacy!($treet);
        end
    end
end
