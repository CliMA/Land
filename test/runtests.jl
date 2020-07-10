using BenchmarkTools
using PlantHydraulics
using Test

PH = PlantHydraulics




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




# FT and NaN test
@testset "Hydraulics --- FT consistency and not NaN" begin
    for FT in [Float32, Float64]
        leaf  = PH.LeafHydraulics{FT}();
        root  = PH.RootHydraulics{FT}();
        stem  = PH.StemHydraulics{FT}();
        roots = [deepcopy(root) for i in 1:5];
        grass = PH.create_grass_like_hs(FT(-2.1), FT(0.5), FT[0,-1,-2,-3], collect(FT,0:1:20));
        palm  = PH.create_palm_like_hs(FT(-2.1), FT(5.5), FT(6), FT[0,-1,-2,-3], collect(FT,0:1:20));
        tree  = PH.create_tree_like_hs(FT(-2.1), FT(5.5), FT(6), FT[0,-1,-2,-3], collect(FT,0:1:20));
        treet = PH.TreeSimple{FT}();
        _vc1 = PH.WeibullSingle{FT}();
        _vc2 = PH.WeibullDual{FT}();
        _sh1 = PH.BrooksCorey{FT}();
        _sh2 = PH.VanGenuchten{FT}();

        # Test the struct
        for data_set in [ leaf, root, stem, roots, grass, palm, tree, treet, _vc1, _vc2, _sh1, _sh2 ]
            recursive_FT_test(data_set, FT)
            recursive_NaN_test(data_set)
        end

        # test xylem_k_ratio function
        _p1  = FT(-2.1)
        _p2  = FT( 0.1)
        _v   = FT(0.9)
        for result in [ PH.xylem_k_ratio(_vc1, _p1),
                        PH.xylem_k_ratio(_vc1, _p2),
                        PH.xylem_k_ratio(_vc1, _p1, _v),
                        PH.xylem_k_ratio(_vc1, _p2, _v),
                        PH.xylem_k_ratio(_vc2, _p1),
                        PH.xylem_k_ratio(_vc2, _p2),
                        PH.xylem_k_ratio(_vc2, _p1, _v),
                        PH.xylem_k_ratio(_vc2, _p2, _v) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        # test the soil functions
        _rwc1 = FT(0.8);
        _rwc2 = FT(0.5);
        for result in [ PH.soil_rwc(_sh1, _p1),
                        PH.soil_rwc(_sh1, _p2),
                        PH.soil_rwc(_sh2, _p1),
                        PH.soil_rwc(_sh2, _p2),
                        PH.soil_k_ratio_rwc(_sh1, _rwc1),
                        PH.soil_k_ratio_rwc(_sh1, _rwc2),
                        PH.soil_k_ratio_rwc(_sh2, _rwc1),
                        PH.soil_k_ratio_rwc(_sh2, _rwc2),
                        PH.soil_k_ratio_p25(_sh1, _p1),
                        PH.soil_p_25(_sh1, _rwc1),
                        PH.soil_p_25(_sh1, _rwc2),
                        PH.soil_p_25(_sh2, _rwc1),
                        PH.soil_p_25(_sh2, _rwc2) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        # test if the soil_p_25 and soil_rwc converge
        _p_giv = FT(-2);
        _r_m_1 = PH.soil_rwc(_sh1, _p_giv);
        _r_m_2 = PH.soil_rwc(_sh2, _p_giv);
        _p_m_1 = PH.soil_p_25(_sh1, _r_m_1);
        _p_m_2 = PH.soil_p_25(_sh2, _r_m_2);
        @test _p_giv ≈ _p_m_1 ≈ _p_m_2;

        # test xylem_p_crit function
        for result in [ PH.xylem_p_crit(_vc1),
                        PH.xylem_p_crit(_vc2) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        # test the xylem_p_from_flow function
        _f_1 = FT(0.001)
        _f_2 = FT(1)
        for result in [ PH.xylem_p_from_flow(leaf, _f_1),
                        PH.xylem_p_from_flow(leaf, _f_2),
                        PH.xylem_p_from_flow(root, _f_1),
                        PH.xylem_p_from_flow(root, _f_2),
                        PH.xylem_p_from_flow(stem, _f_1),
                        PH.xylem_p_from_flow(stem, _f_2) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        # test the leaf_xylem_risk function
        for result in [ PH.leaf_xylem_risk(leaf, _f_1),
                        PH.leaf_xylem_risk(leaf, _f_2) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        # test the root q functions
        _p1 = FT(0)
        _p2 = FT(-0.5)
        _p3 = FT(-1.0)
        for result in [ PH.root_q_from_pressure(root, _p1),
                        PH.root_q_from_pressure(root, _p2),
                        PH.root_q_from_pressure(root, _p3) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        # test the root q functions
        _p1 = FT(0)
        _p2 = FT(-0.5)
        _p3 = FT(-1.0)
        for result in [ PH.root_qs_p_from_q(roots, FT(0)),
                        PH.root_qs_p_from_q(roots, _f_1),
                        PH.root_qs_p_from_q(roots, _f_2) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        # test the leaf_xylem_risk function
        for result in [ PH.leaf_e_crit(leaf) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        # test the xylem_p_from_flow function for plants
        for result in [ PH.xylem_p_from_flow(grass, _f_1),
                        PH.xylem_p_from_flow(grass, _f_2),
                        PH.xylem_p_from_flow(palm, _f_1),
                        PH.xylem_p_from_flow(palm, _f_2),
                        PH.xylem_p_from_flow(tree, _f_1),
                        PH.xylem_p_from_flow(tree, _f_2),
                        PH.xylem_p_from_flow(treet, _f_1),
                        PH.xylem_p_from_flow(treet, _f_2) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        # test the tree_e_crit function
        for result in [ PH.tree_e_crit(grass),
                        PH.tree_e_crit(palm),
                        PH.tree_e_crit(tree),
                        PH.tree_e_crit(treet) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end
    end
end
