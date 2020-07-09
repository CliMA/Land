# include("./scripts/leaf investment/yujie_111_map_sun_shade.jl")

node   = Yujie111Init();
Yujie111UpdateLeaf(node, 750.0, 100.0)
Yujie111UpdateSoilFromP(node, 0.0);
envir  = [25.0 40.0 1.5 2.0];
zenith = 30.0
r_all  = 1000.0
matrix = Yujie111GetOptimalFsMap(node, envir, zenith, r_all)

writedlm("./sun_shade.txt", matrix)