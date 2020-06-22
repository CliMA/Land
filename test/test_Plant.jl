# FT and NaN tests
@testset "Plant --- FT consistency and not NaN" begin
    for FT in [Float32, Float64]
        for data_set in [ PT.ESMBallBerry{FT}(),
                          PT.ESMGentine{FT}(),
                          PT.ESMLeuning{FT}(),
                          PT.ESMMedlyn{FT}(),
                          PT.OSMWAP{FT}(),
                          PT.OSMWAPMod{FT}(),
                          PT.RootLayer{FT}(),
                          PT.create_root_list(FT,4,-1.0),
                          PT.Stem{FT}(),
                          PT.create_branch_list(FT,20,4.0,8.0),
                          PT.Leaf{FT}(),
                          PT.CanopyLayer{FT,325}(),
                          PT.Tree{FT,5,20,325}() ]
            recursive_FT_test(data_set, FT)
            recursive_NaN_test(data_set, FT)
        end
    end
end




# Hydraulic and stomatal functions
@testset "Plant --- Hydraulic and stomatal functions" begin
    for FT in [Float32, Float64]
        tree = PT.Tree{FT,5,20,325}()
        PT.update_tree_e_crit!(tree)
        PT.update_leaf_ak_max!(tree)
        PT.update_struct_from_q!(tree.trunk, FT(0.01)),
        PT.update_struct_from_q!(tree.branch_list[1], FT(0.01)),
        PT.update_struct_from_q!(tree.canopy_list[1].leaf_list[1], FT(0.0001))
        PT.update_tree_with_time!(tree, FT(1.0), PT.OSMEller())
        PT.update_tree_with_time!(tree, FT(1.0), PT.OSMSperry())
        PT.update_tree_with_time!(tree, FT(1.0), PT.OSMWang())
        PT.update_tree_with_time!(tree, FT(1.0), PT.OSMWAP{FT}())
        PT.update_tree_with_time!(tree, FT(1.0), PT.OSMWAPMod{FT}())
        PT.update_tree_with_time!(tree, FT(1.0), PT.ESMBallBerry{FT}())
        PT.update_tree_with_time!(tree, FT(1.0), PT.ESMGentine{FT}())
        PT.update_tree_with_time!(tree, FT(1.0), PT.ESMLeuning{FT}())
        PT.update_tree_with_time!(tree, FT(1.0), PT.ESMMedlyn{FT}())
        for result in [ PT.get_p_base_q_list_from_q(tree, FT(1.0)),
                        PT.get_q_layer_from_p_base(tree.root_list[1], FT(-0.5)),
                        PT.get_struct_p_end_from_q(tree.trunk, FT(0.1)),
                        PT.get_struct_p_end_from_q(tree.branch_list[1], FT(0.1)),
                        PT.get_struct_p_end_from_q(tree.canopy_list[1].leaf_list[1], FT(0.001)),
                        PT.get_marginal_gain(tree.canopy_list[1], 1, tree.photo_para_set),
                        PT.get_marginal_gain(tree.canopy_list[1], tree.photo_para_set),
                        PT.get_marginal_penalty(tree.canopy_list[1], 1, PT.OSMEller()),
                        PT.get_marginal_penalty(tree.canopy_list[1], PT.OSMEller()),
                        PT.get_marginal_penalty(tree.canopy_list[1], 1, PT.OSMSperry()),
                        PT.get_marginal_penalty(tree.canopy_list[1], PT.OSMSperry()),
                        PT.get_marginal_penalty(tree.canopy_list[1], 1, PT.OSMWang()),
                        PT.get_marginal_penalty(tree.canopy_list[1], PT.OSMWang()),
                        PT.get_marginal_penalty(tree.canopy_list[1], 1, PT.OSMWAP{FT}()),
                        PT.get_marginal_penalty(tree.canopy_list[1], PT.OSMWAP{FT}()),
                        PT.get_marginal_penalty(tree.canopy_list[1], 1, PT.OSMWAPMod{FT}()),
                        PT.get_marginal_penalty(tree.canopy_list[1], PT.OSMWAPMod{FT}()),
                        PT.get_empirical_gsw_from_model(PT.ESMBallBerry{FT}(), FT(5.0), FT(101325.0), FT(40.0), FT(0.5), FT(1500), FT(4.5), FT(0.5)),
                        PT.get_empirical_gsw_from_model(PT.ESMGentine{FT}(), FT(5.0), FT(101325.0), FT(40.0), FT(0.5)),
                        PT.get_empirical_gsw_from_model(PT.ESMLeuning{FT}(), FT(5.0), FT(101325.0), FT(40.0), FT(0.5), FT(1500), FT(4.5), FT(0.5)),
                        PT.get_empirical_gsw_from_model(PT.ESMMedlyn{FT}(), FT(5.0), FT(101325.0), FT(40.0), FT(0.5), FT(1500), FT(4.5), FT(0.5)),
                        PT.get_empirical_gsw_from_model(PT.ESMBallBerry{FT}(), ones(FT,10).*5, FT(101325.0), FT(40.0), FT(0.5), ones(FT,10)*1500, ones(FT,10).*4, ones(FT,10)/2),
                        PT.get_empirical_gsw_from_model(PT.ESMGentine{FT}(), ones(FT,10).*5, FT(101325.0), FT(40.0), ones(FT,10)/2),
                        PT.get_empirical_gsw_from_model(PT.ESMLeuning{FT}(), ones(FT,10).*5, FT(101325.0), FT(40.0), FT(0.5), ones(FT,10)*1500, ones(FT,10).*4, ones(FT,10)/2),
                        PT.get_empirical_gsw_from_model(PT.ESMMedlyn{FT}(), ones(FT,10).*5, FT(101325.0), FT(40.0), FT(0.5), ones(FT,10)*1500, ones(FT,10).*4, ones(FT,10)/2),
                         ]
            recursive_FT_test(result, FT)
            recursive_NaN_test(result, FT)
        end 
    end
end
