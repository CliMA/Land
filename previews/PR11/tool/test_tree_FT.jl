# include("docs/src/tool/test_tree_FT.jl")

using Land.Plant

FT = Plant.FT

tree = Tree{FT}();

p_base,q_list = Plant.get_p_base_q_list_from_q(tree, FT(1.0))
println(typeof(p_base))
println(typeof(q_list))

q_layer = Plant.get_q_layer_from_p_base(tree.roots.root_list[1], FT(-0.1));
println(typeof(q_layer))

p_end = Plant.get_struct_p_end_from_q(tree.roots.root_list[1], FT(0.1); p_ini=FT(Inf));
println(typeof(p_end))

p_end = Plant.get_struct_p_end_from_q(tree.trunk, FT(1.0); p_ini=FT(Inf));
println(typeof(p_end))

p_end = Plant.get_struct_p_end_from_q(tree.branch.branch_list[1], FT(0.1); p_ini=FT(Inf));
println(typeof(p_end))

p_end = Plant.get_struct_p_end_from_q(tree.canopy.canopy_list[1].leaf_list[1], FT(0.05); p_ini=FT(Inf));
println(typeof(p_end))

temp = Plant.get_leaf_an_ag_r_from_pi();
println(typeof(temp))

temp = Plant.get_leaf_an_ag_r_pi_from_gsc();
println(typeof(temp))

temp = Plant.get_leaf_an_ag_r_pi_from_gsc_list();
println(typeof(temp))

temp = Plant.get_leaf_j(FT(100.0), FT(1000.0));
println(typeof(temp))

temp = Plant.get_leaf_jmax(FT(100.0), FT(298.0));
println(typeof(temp))

temp = Plant.get_leaf_r_from_r25(FT(2.0), FT(298.0));
println(typeof(temp))

temp = Plant.get_leaf_r_from_v25(FT(100.0), FT(298.0));
println(typeof(temp))

temp = Plant.get_leaf_vcmax(FT(100.0), FT(298.0));
println(typeof(temp))

temp = Plant.get_marginal_gain(tree.canopy.canopy_list[1], 1);
println(typeof(temp))

temp = Plant.get_marginal_gain(tree.canopy.canopy_list[1]);
println(typeof(temp))

temp = Plant.get_marginal_penalty_wang(tree.canopy.canopy_list[1], 1);
println(typeof(temp))

temp = Plant.get_marginal_penalty_wang(tree.canopy.canopy_list[1]);
println(typeof(temp))

temp = Plant.get_a_par_curve();
println(typeof(temp))

temp = Plant.get_a_pi_curve();
println(typeof(temp))

temp = Plant.get_relative_surface_tension(FT(298.0));
println(typeof(temp))

temp = Plant.get_relative_viscosity(FT(298.0));
println(typeof(temp))

temp = Plant.get_saturated_vapor_pressure(FT(298.0));
println(typeof(temp))

temp = Plant.get_saturated_vapor_pressure(ones(FT,10) .* FT(298.0));
println(typeof(temp))

temp = Plant.get_specific_latent_heat(FT(298.5));
println(typeof(temp))
