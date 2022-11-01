# This file benchmarks the functions one by one
using BenchmarkTools
using ClimaCache
using Photosynthesis


FT = Float64;


# File colimit
for met in [ClimaCache.MinimumColimit{FT}(),
            ClimaCache.QuadraticColimit{FT}(),
            ClimaCache.SerialColimit{FT}(),
            ClimaCache.SquareColimit{FT}()]
    val1 = FT(1);
    val2 = FT(2);
    @btime Photosynthesis.colimited_rate($val1, $val2, $met);
end

for psm in [ClimaCache.C3CytochromeModel{FT}(),
            ClimaCache.C3VJPModel{FT}(),
            ClimaCache.C4VJPModel{FT}()]
    beta = FT(0.5);
    @btime Photosynthesis.colimit_photosynthesis!($psm);
    @btime Photosynthesis.colimit_photosynthesis!($psm; β = $beta);
end


# File etr
begin
    psm1 = ClimaCache.C3CytochromeModel{FT}();
    psm2 = ClimaCache.C3VJPModel{FT}();
    psm3 = ClimaCache.C4VJPModel{FT}();
    prc1 = ClimaCache.CytochromeReactionCenter{FT}();
    prc2 = ClimaCache.VJPReactionCenter{FT}();
    prc3 = ClimaCache.VJPReactionCenter{FT}();
    ppar = FT(100);
    p_i  = FT(20);
    beta = FT(0.5);
    @btime Photosynthesis.photosystem_electron_transport!($psm1, $prc1, $ppar, $p_i);
    @btime Photosynthesis.photosystem_electron_transport!($psm2, $prc2, $ppar, $p_i);
    @btime Photosynthesis.photosystem_electron_transport!($psm3, $prc3, $ppar, $p_i);
    @btime Photosynthesis.photosystem_electron_transport!($psm1, $prc1, $ppar, $p_i; β = $beta);
    @btime Photosynthesis.photosystem_electron_transport!($psm2, $prc2, $ppar, $p_i; β = $beta);
    @btime Photosynthesis.photosystem_electron_transport!($psm3, $prc3, $ppar, $p_i; β = $beta);
end;


# File fluorescence
begin
    psm1 = ClimaCache.C3CytochromeModel{FT}();
    psm2 = ClimaCache.C3VJPModel{FT}();
    psm3 = ClimaCache.C4VJPModel{FT}();
    prc1 = ClimaCache.CytochromeReactionCenter{FT}();
    prc2 = ClimaCache.VJPReactionCenter{FT}();
    prc3 = ClimaCache.VJPReactionCenter{FT}();
    ppar = FT(100);
    beta = FT(0.5);
    @btime Photosynthesis.photosystem_coefficients!($psm1, $prc1, $ppar);
    @btime Photosynthesis.photosystem_coefficients!($psm2, $prc2, $ppar);
    @btime Photosynthesis.photosystem_coefficients!($psm3, $prc3, $ppar);
    @btime Photosynthesis.photosystem_coefficients!($psm1, $prc1, $ppar; β = $beta);
    @btime Photosynthesis.photosystem_coefficients!($psm2, $prc2, $ppar; β = $beta);
    @btime Photosynthesis.photosystem_coefficients!($psm3, $prc3, $ppar; β = $beta);
end;


# File light_limited
begin
    psm1 = ClimaCache.C3CytochromeModel{FT}();
    psm2 = ClimaCache.C3VJPModel{FT}();
    psm3 = ClimaCache.C4VJPModel{FT}();
    prc1 = ClimaCache.CytochromeReactionCenter{FT}();
    prc2 = ClimaCache.VJPReactionCenter{FT}();
    prc3 = ClimaCache.VJPReactionCenter{FT}();
    air  = ClimaCache.AirLayer{FT}();
    g_lc = FT(0.1);
    beta = FT(0.5);
    @btime Photosynthesis.light_limited_rate!($psm1);
    @btime Photosynthesis.light_limited_rate!($psm2);
    @btime Photosynthesis.light_limited_rate!($psm3);
    @btime Photosynthesis.light_limited_rate!($psm1, $prc1, $air, $g_lc; β = $beta);
    @btime Photosynthesis.light_limited_rate!($psm2, $prc2, $air, $g_lc; β = $beta);
    @btime Photosynthesis.light_limited_rate!($psm3, $prc3, $air, $g_lc; β = $beta);
end;


# File product_limited
begin
    psm1 = ClimaCache.C3CytochromeModel{FT}();
    psm2 = ClimaCache.C3VJPModel{FT}();
    psm3 = ClimaCache.C4VJPModel{FT}();
    prc1 = ClimaCache.CytochromeReactionCenter{FT}();
    prc2 = ClimaCache.VJPReactionCenter{FT}();
    prc3 = ClimaCache.VJPReactionCenter{FT}();
    air  = ClimaCache.AirLayer{FT}();
    g_lc = FT(0.1);
    p_i  = FT(20);
    beta = FT(0.5);
    @btime Photosynthesis.product_limited_rate!($psm1, $p_i);
    @btime Photosynthesis.product_limited_rate!($psm2, $p_i);
    @btime Photosynthesis.product_limited_rate!($psm3, $p_i);
    @btime Photosynthesis.product_limited_rate!($psm1, $p_i; β = $beta);
    @btime Photosynthesis.product_limited_rate!($psm2, $p_i; β = $beta);
    @btime Photosynthesis.product_limited_rate!($psm3, $p_i; β = $beta);
    @btime Photosynthesis.product_limited_rate!($psm1, $air, $g_lc);
    @btime Photosynthesis.product_limited_rate!($psm2, $air, $g_lc);
    @btime Photosynthesis.product_limited_rate!($psm3, $air, $g_lc);
    @btime Photosynthesis.product_limited_rate!($psm1, $air, $g_lc; β = $beta);
    @btime Photosynthesis.product_limited_rate!($psm2, $air, $g_lc; β = $beta);
    @btime Photosynthesis.product_limited_rate!($psm3, $air, $g_lc; β = $beta);
end;


# File rubisco_limited
begin
    psm1 = ClimaCache.C3CytochromeModel{FT}();
    psm2 = ClimaCache.C3VJPModel{FT}();
    psm3 = ClimaCache.C4VJPModel{FT}();
    prc1 = ClimaCache.CytochromeReactionCenter{FT}();
    prc2 = ClimaCache.VJPReactionCenter{FT}();
    prc3 = ClimaCache.VJPReactionCenter{FT}();
    air  = ClimaCache.AirLayer{FT}();
    g_lc = FT(0.1);
    p_i  = FT(20);
    beta = FT(0.5);
    @btime Photosynthesis.rubisco_limited_rate!($psm1, $p_i);
    @btime Photosynthesis.rubisco_limited_rate!($psm2, $p_i);
    @btime Photosynthesis.rubisco_limited_rate!($psm3, $p_i);
    @btime Photosynthesis.rubisco_limited_rate!($psm1, $p_i; β = $beta);
    @btime Photosynthesis.rubisco_limited_rate!($psm2, $p_i; β = $beta);
    @btime Photosynthesis.rubisco_limited_rate!($psm3, $p_i; β = $beta);
    @btime Photosynthesis.rubisco_limited_rate!($psm1, $air, $g_lc);
    @btime Photosynthesis.rubisco_limited_rate!($psm2, $air, $g_lc);
    @btime Photosynthesis.rubisco_limited_rate!($psm3, $air, $g_lc);
    @btime Photosynthesis.rubisco_limited_rate!($psm1, $air, $g_lc; β = $beta);
    @btime Photosynthesis.rubisco_limited_rate!($psm2, $air, $g_lc; β = $beta);
    @btime Photosynthesis.rubisco_limited_rate!($psm3, $air, $g_lc; β = $beta);
end;


# temperature file
for tpd in [ClimaCache.Arrhenius{FT}(T_REF = 298, VAL_REF = 40, ΔHA = 80000),
            ClimaCache.ArrheniusPeak{FT}(T_REF = 298, VAL_REF = 40 , ΔHA = 50000, ΔHD = 400000, ΔSV = 1000),
            ClimaCache.Q10{FT}(Q_10 = 1.4, T_REF = 298, VAL_REF = 1),
            ClimaCache.Q10Peak{FT}(Q_10 = 1.4, T_REF = 298, VAL_REF = 1, ΔHD = 400000, ΔSV = 1000)]
    tem = FT(300);
    ref = FT(310);
    @btime Photosynthesis.temperature_correction($tpd, $tem);
    @btime Photosynthesis.temperature_correction($tpd, $tem; t_ref = $ref);
    @btime Photosynthesis.temperature_corrected_value($tpd, $tem);
    @btime Photosynthesis.temperature_corrected_value($tpd, $tem; t_ref = $ref);
    @btime Photosynthesis.∂R∂T($tpd, $ref, $tem);
end
for psm in [ClimaCache.C3CytochromeModel{FT}(),
            ClimaCache.C3VJPModel{FT}(),
            ClimaCache.C4VJPModel{FT}()]
    air  = ClimaCache.AirLayer{FT}();
    tem = FT(300);
    @btime Photosynthesis.rubisco_limited_rate!($psm, $air, $tem);
    @btime Photosynthesis.∂R∂T($psm, $tem);
end;
for lfv in [ClimaCache.Leaf{FT}(),
            ClimaCache.Leaves1D{FT}(),
            ClimaCache.Leaves2D{FT}()]
    @btime Photosynthesis.∂R∂T($lfv);
end


# model file
for lfv in [ClimaCache.Leaf{FT}(),
            ClimaCache.Leaves1D{FT}(),
            ClimaCache.Leaves2D{FT}()]
    air  = ClimaCache.AirLayer{FT}();
    modg = ClimaCache.GCO₂Mode();
    modp = ClimaCache.PCO₂Mode();
    beta = ClimaCache.BetaFunction{FT}();
    betg = ClimaCache.BetaParameterG1();
    betv = ClimaCache.BetaParameterVcmax();
    g_lc = FT(0.1);
    ppar = FT(100);
    tem  = FT(300);
    corr = FT(0.5);
    @btime Photosynthesis.leaf_photosynthesis!($lfv, $air, $g_lc, $ppar, $tem);
    for mode in [ClimaCache.GCO₂Mode(), ClimaCache.PCO₂Mode()]
        @btime Photosynthesis.leaf_photosynthesis!($lfv, $air, $mode);
        @btime Photosynthesis.leaf_photosynthesis!($lfv, $air, $mode, $corr);
        # for stm in [ClimaCache.AndereggSM{FT}(),
        #             ClimaCache.EllerSM{FT}(),
        #             ClimaCache.SperrySM{FT}(),
        #             ClimaCache.WangSM{FT}(),
        #             ClimaCache.Wang2SM{FT}(),
        #             ClimaCache.BallBerrySM{FT}(),
        #             ClimaCache.GentineSM{FT}(),
        #             ClimaCache.LeuningSM{FT}(),
        #             ClimaCache.MedlynSM{FT}()]
        #     @btime Photosynthesis.leaf_photosynthesis!($lfv, $air, $mode, $stm);
        # end;
        @btime Photosynthesis.leaf_photosynthesis!($lfv, $air, $mode, $beta, $betg);
        @btime Photosynthesis.leaf_photosynthesis!($lfv, $air, $mode, $beta, $betv);
    end;
end
for spc in [ClimaCache.MonoElementSPAC{FT}(),
            ClimaCache.MonoMLGrassSPAC{FT}(),
            ClimaCache.MonoMLPalmSPAC{FT}(),
            ClimaCache.MonoMLTreeSPAC{FT}()]
    for mode in [ClimaCache.GCO₂Mode(), ClimaCache.PCO₂Mode()]
        @btime Photosynthesis.leaf_photosynthesis!($spc, $mode);
    end;
end;
