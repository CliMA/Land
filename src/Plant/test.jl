using BenchmarkTools




include("include_plant.jl");
tree = StructTree();
@time update_tree_e_crit(tree);




@btime update_tree_with_time(tree,1.0,updating=false);




visualize_struct_tree(tree)



