# include("scripts/yujie_111_price_of_water_1.jl")

# pre-drought
node   = Yujie111Init();
Yujie111UpdateLeaf(node, 750.0, 100.0)
Yujie111UpdateSoilFromP(node, 0.0);
envir  = [25.0 40.0 1.5 2.0];
zenith = 30.0
r_all  = 1000.0
matrix = Yujie111EvaluateModel(node, envir, zenith, r_all)
writedlm("./price_of_water_1.txt", matrix)




# drought
node   = Yujie111Init();
Yujie111UpdateLeaf(node, 750.0, 100.0)
Yujie111UpdateSoilFromP(node, 2.2);
Yujie111UpdateLegacy(node, 0.0, 0.0, 0.5);
envir  = [25.0 40.0 1.5 2.0];
zenith = 30.0
r_all  = 1000.0
matrix = Yujie111EvaluateModel(node, envir, zenith, r_all)
writedlm("./price_of_water_2.txt", matrix)




# post-drought
node   = Yujie111Init();
Yujie111UpdateLeaf(node, 750.0, 100.0)
Yujie111UpdateSoilFromP(node, 2.2);
Yujie111UpdateLegacy(node, 0.0, 0.0, 0.5);
Yujie111UpdateSoilFromP(node, 0.0);
envir  = [25.0 40.0 1.5 2.0];
zenith = 30.0
r_all  = 1000.0
matrix = Yujie111EvaluateModel(node, envir, zenith, r_all)
writedlm("./price_of_water_3.txt", matrix)
