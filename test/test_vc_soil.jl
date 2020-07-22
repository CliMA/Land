# FT and NaN test
println("\nTesting and benchmarking the soil VC functions...")
@testset "Hydraulics --- soil VC" begin
    for FT in [Float32, Float64]
        _sh1  = BrooksCorey{FT}();
        _sh2  = VanGenuchten{FT}();

        # test soil functions
        _rwc = FT(0.5);
        _p   = FT(-0.5);
        for result in [ soil_rwc(_sh1, _p),
                        soil_rwc(_sh2, _p),
                        soil_k_ratio_rwc(_sh1, _rwc),
                        soil_k_ratio_rwc(_sh2, _rwc),
                        soil_k_ratio_swc(_sh1, _rwc*_sh1.Θs),
                        soil_k_ratio_swc(_sh2, _rwc*_sh2.Θs),
                        soil_k_ratio_p25(_sh1, _p),
                        soil_k_ratio_p25(_sh2, _p),
                        soil_p_25_rwc(_sh1, _rwc),
                        soil_p_25_rwc(_sh2, _rwc),
                        soil_p_25_swc(_sh1, _rwc*_sh1.Θs),
                        soil_p_25_swc(_sh2, _rwc*_sh2.Θs) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        # test if the soil_p_25 and soil_rwc converge
        _p_giv = FT(-2);
        _r_m_1 = soil_erwc(_sh1, _p_giv);
        _r_m_2 = soil_erwc(_sh2, _p_giv);
        _p_m_1 = soil_p_25_erwc(_sh1, _r_m_1);
        _p_m_2 = soil_p_25_erwc(_sh2, _r_m_2);
        @test _p_giv ≈ _p_m_1 ≈ _p_m_2;

        _r_m_1 = soil_rwc(_sh1, _p_giv);
        _r_m_2 = soil_rwc(_sh2, _p_giv);
        _p_m_1 = soil_p_25_rwc(_sh1, _r_m_1);
        _p_m_2 = soil_p_25_rwc(_sh2, _r_m_2);
        @test _p_giv ≈ _p_m_1 ≈ _p_m_2;

        _r_m_1 = soil_swc(_sh1, _p_giv);
        _r_m_2 = soil_swc(_sh2, _p_giv);
        _p_m_1 = soil_p_25_swc(_sh1, _r_m_1);
        _p_m_2 = soil_p_25_swc(_sh2, _r_m_2);
        @test _p_giv ≈ _p_m_1 ≈ _p_m_2;

        if benchmarking
            @btime soil_rwc($_sh1, $_p);
            @btime soil_rwc($_sh1, $_p);
            @btime soil_k_ratio_rwc($_sh1, $_rwc);
            @btime soil_k_ratio_rwc($_sh2, $_rwc);
            @btime soil_k_ratio_swc($_sh1, $_rwc * $_sh1.Θs);
            @btime soil_k_ratio_swc($_sh2, $_rwc * $_sh2.Θs);
            @btime soil_k_ratio_p25($_sh1, $_p);
            @btime soil_k_ratio_p25($_sh2, $_p);
            @btime soil_p_25_rwc($_sh1, $_rwc);
            @btime soil_p_25_rwc($_sh2, $_rwc);
            @btime soil_p_25_swc($_sh1, $_rwc * $_sh1.Θs);
            @btime soil_p_25_swc($_sh2, $_rwc * $_sh2.Θs);
        end
    end
end
