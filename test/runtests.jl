using BenchmarkTools
using CSV
using DataFrames
using Plants
using Photosynthesis
using Test




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




# test the structs
println("\nTesting the FT and NaN of the structs...")
@testset "FT and NaN --- Types" begin
    for FT in [Float32, Float64]
        node   = SPACSimple{FT}();
        cont1L = SPACContainer1L{FT}();
        cont2L = SPACContainer2L{FT}();

        for data in [node, cont1L, cont2L]
            recursive_FT_test(data, FT);
            recursive_NaN_test(data);
        end
    end
end




# test and benchmark the big_leaf_partition!
println("\nTesting and Benchmarking the big_leaf_partition!...")
@testset "Testing + Benchmarking --- big_leaf_partition!" begin
    for FT in [Float32, Float64]
        node   = SPACSimple{FT}();
        zenith = FT(30);
        r_all  = FT(1000);

        big_leaf_partition!(node, zenith, r_all);
        recursive_FT_test(node, FT);
        recursive_NaN_test(node);

        @btime big_leaf_partition!($node, $zenith, $r_all);
    end
end




# test the gain_risk_map
println("\nTesting the gain_risk_map Function...")
@testset "Testing --- gain_risk_map" begin
    for FT in [Float32, Float64]
        node   = SPACSimple{FT}();
        photo  = C3CLM(FT);
        zenith = FT(30);
        r_all  = FT(1000);

        big_leaf_partition!(node, zenith, r_all);
        mat = gain_risk_map(node, photo);
        recursive_FT_test(mat, FT);
        recursive_NaN_test(mat);
    end
end




# test and benchmark the leaf_gas_exchange_nonopt!
println("\nTesting and Benchmarking the leaf_gas_exchange_nonopt! Functions...")
@testset "Testing + Benchmarking --- leaf_gas_exchange_nonopt!" begin
    for FT in [Float32, Float64]
        node   = SPACSimple{FT}();
        photo  = C3CLM(FT);
        zenith = FT(30);
        r_all  = FT(1000);
        flow   = FT(4);
        f_sl   = FT(2.5);
        f_sh   = FT(1.5);

        big_leaf_partition!(node, zenith, r_all);
        leaf_gas_exchange_nonopt!(node, photo, flow);
        recursive_FT_test(node, FT);
        recursive_NaN_test(node);
        leaf_gas_exchange_nonopt!(node, photo, f_sl, f_sh);
        recursive_FT_test(node, FT);
        recursive_NaN_test(node);

        @btime leaf_gas_exchange_nonopt!($node, $photo, $flow);
        @btime leaf_gas_exchange_nonopt!($node, $photo, $f_sl, $f_sh);
    end
end




# test and benchmark the leaf_gas_exchange!
println("\nTesting and Benchmarking the leaf_gas_exchange! Functions...")
@testset "Testing + Benchmarking --- leaf_gas_exchange!" begin
    for FT in [Float32, Float64]
        node   = SPACSimple{FT}();
        photo  = C3CLM(FT);
        zenith = FT(30);
        r_all  = FT(1000);
        flow   = FT(4);
        f_sl   = FT(2.5);
        f_sh   = FT(1.5);

        big_leaf_partition!(node, zenith, r_all);
        leaf_gas_exchange!(node, photo, flow);
        recursive_FT_test(node, FT);
        recursive_NaN_test(node);
        leaf_gas_exchange!(node, photo, f_sl, f_sh);
        recursive_FT_test(node, FT);
        recursive_NaN_test(node);

        @btime leaf_gas_exchange!($node, $photo, $flow);
        @btime leaf_gas_exchange!($node, $photo, $f_sl, $f_sh);
    end
end




# test and benchmark the leaf_temperature*
println("\nTesting and Benchmarking the leaf_temperature* Functions...")
@testset "Testing + Benchmarking --- leaf_temperature*" begin
    for FT in [Float32, Float64]
        node = SPACSimple{FT}();
        rad  = FT(300);
        flow = FT(4);

        for result in [ leaf_temperature(node, rad, flow),
                        leaf_temperature_shaded(node, rad, flow),
                        leaf_temperature_sunlit(node, rad, flow) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        @btime leaf_temperature($node, $rad, $flow);
        @btime leaf_temperature_shaded($node, $rad, $flow);
        @btime leaf_temperature_sunlit($node, $rad, $flow);
    end
end




# test and benchmark the optimize_flows!
println("\nTesting and Benchmarking the optimize_flows! Functions...")
@testset "Testing + Benchmarking --- optimize_flows!" begin
    for FT in [Float32, Float64]
        node   = SPACSimple{FT}();
        photo  = C3CLM(FT);
        zenith = FT(30);
        r_all  = FT(1000);

        big_leaf_partition!(node, zenith, r_all);
        optimize_flows!(node, photo);
        recursive_FT_test(node, FT);
        recursive_NaN_test(node);

        # reset the node before benchmarking
        node   = SPACSimple{FT}();
        @btime optimize_flows!($node, $photo);
    end
end




# test and benchmark the atmosheric* functions
println("\nTesting and Benchmarking the atmosheric* Functions...")
@testset "Testing + Benchmarking --- atmosheric*" begin
    for FT in [Float32, Float64]
        h = FT(1000);

        for result in [ atmospheric_pressure(h),
                        atmospheric_pressure_ratio(h),
                        ppm_to_Pa(h) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        @btime atmospheric_pressure($h);
        @btime atmospheric_pressure_ratio($h);
        @btime ppm_to_Pa($h);
    end
end




# test and benchmark the zenith_angle
println("\nTesting and Benchmarking the zenith_angle Functions...")
@testset "Testing + Benchmarking --- zenith_angle" begin
    for FT in [Float32, Float64]
        latd = FT(10);
        decd = FT(10);
        lhad = FT(10);
        day  = FT(100)
        hour = FT(13)
        minu = FT(30)

        for result in [ zenith_angle(latd, decd, lhad),
                        zenith_angle(latd, day, hour),
                        zenith_angle(latd, day, hour, minu) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        @btime zenith_angle($latd, $decd, $lhad);
        @btime zenith_angle($latd, $day, $hour);
        @btime zenith_angle($latd, $day, $hour, $minu);
    end
end




# test and benchmark the annual_profit
println("\nTesting and Benchmarking the annual_profit Functions...")
@testset "Testing + Benchmarking --- annual_profit" begin
    weat = DataFrame!(CSV.File("../data/gs_sample.csv"));
    for FT in [Float32, Float64]
        node    = SPACSimple{FT}();
        photo   = C3CLM(FT);
        weatmat = Matrix{FT}(weat);

        gscp = annual_profit(node, photo, weatmat);
        recursive_FT_test(gscp, FT);
        recursive_NaN_test(gscp);

        # reset the node before benchmarking
        node = SPACSimple{FT}();
        @btime annual_profit($node, $photo, $weatmat);
    end
end




# test and benchmark the annual_simulation!
println("\nTesting annual_simulation! Functions...")
@testset "Testing + Benchmarking --- annual_simulation!" begin
    weat = DataFrame!(CSV.File("../data/gs_sample.csv"));
    for FT in [Float32, Float64]
        node  = SPACSimple{FT}();
        photo = C3CLM(FT);
        df    = DataFrame();
        df[!, "Time"  ]  = weat.Day + weat.Hour / FT(24);
        df[!, "T_air" ] .= FT(0);
        df[!, "D_air" ] .= FT(0);
        df[!, "Wind"  ] .= FT(0);
        df[!, "Rain"  ] .= FT(0);
        df[!, "Ca"    ] .= FT(0);
        df[!, "SWC"   ] .= FT(0);
        df[!, "P_soil"] .= FT(0);
        df[!, "H_sun" ] .= FT(0);
        df[!, "A_net" ] .= FT(0);
        df[!, "E_crit"] .= FT(0);
        df[!, "LAI_sl"] .= FT(0);
        df[!, "PAR_sl"] .= FT(0);
        df[!, "RAD_sl"] .= FT(0);
        df[!, "E_sl"  ] .= FT(0);
        df[!, "P_sl"  ] .= FT(0);
        df[!, "An_sl" ] .= FT(0);
        df[!, "Ag_sl" ] .= FT(0);
        df[!, "C_sl"  ] .= FT(0);
        df[!, "G_sl"  ] .= FT(0);
        df[!, "T_sl"  ] .= FT(0);
        df[!, "LAI_sh"] .= FT(0);
        df[!, "PAR_sh"] .= FT(0);
        df[!, "RAD_sh"] .= FT(0);
        df[!, "E_sh"  ] .= FT(0);
        df[!, "P_sh"  ] .= FT(0);
        df[!, "An_sh" ] .= FT(0);
        df[!, "Ag_sh" ] .= FT(0);
        df[!, "C_sh"  ] .= FT(0);
        df[!, "G_sh"  ] .= FT(0);
        df[!, "T_sh"  ] .= FT(0);

        annual_simulation!(node, photo, weat, df);
        recursive_FT_test(node, FT);
        recursive_NaN_test(node);
    end
end




# test and benchmark the leaf_allocation!
println("\nTesting and Benchmarking the leaf_allocation! Functions...")
@testset "Testing + Benchmarking --- leaf_allocation!" begin
    for FT in [Float32, Float64]
        node  = SPACSimple{FT}();
        photo = C3CLM(FT);
        laba  = FT(1000);
        vmax  = FT(80);

        leaf_allocation!(node, laba);
        recursive_FT_test(node, FT);
        recursive_NaN_test(node);
        leaf_allocation!(node, photo, vmax);
        recursive_FT_test(node, FT);
        recursive_NaN_test(node);
        leaf_allocation!(node, photo, laba, vmax);
        recursive_FT_test(node, FT);
        recursive_NaN_test(node);

        @btime leaf_allocation!($node, $laba);
        @btime leaf_allocation!($node, $photo, $vmax);
        @btime leaf_allocation!($node, $photo, $laba, $vmax);
    end
end




# test and benchmark the optimize_leaf
println("\nTesting and Benchmarking the optimize_leaf! Functions...")
@testset "Testing + Benchmarking --- optimize_leaf!" begin
    weat = DataFrame!(CSV.File("../data/gs_sample.csv"));
    for FT in [Float32, Float64]
        node    = SPACSimple{FT}();
        photo   = C3CLM(FT);
        weatmat = Matrix{FT}(weat);

        optimize_leaf!(node, photo, weatmat; displaying=false);
        recursive_FT_test(node, FT);
        recursive_NaN_test(node);

        # reset the node before benchmarking
        node = SPACSimple{FT}();
        @btime optimize_leaf!($node, $photo, $weatmat; displaying=false);
    end
end
