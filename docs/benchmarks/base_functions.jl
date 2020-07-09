using BenchmarkTools
using Photosynthesis
using Plants

FT = Float32




#=
# Benchmarking functions in the Earth Folder
println("Benchmarking the Earth sub-module...");
height = FT(1000);
latd   = FT(30);
decd   = FT(10);
lhad   = FT(15);
minute = FT(10);
@btime atmospheric_pressure_ratio($height);
@btime atmospheric_pressure($height);
@btime ppm_to_Pa($height);
@btime zenith_angle($latd, $decd, $lhad);
@btime zenith_angle($latd, 100, 11);
@btime zenith_angle($latd, 100, 11, $minute);




# Benchmarking functions in the Radiation Folder
println("Benchmarking the Radiation sub-module...");
partition = zeros(FT,6);
zenith    = FT(30);
lai       = FT(3);
r_all     = FT(1000);
@btime big_leaf_partition!($partition, $zenith, $lai, $r_all);




# Benchmarking functions in the Tree Folder
println("Benchmarking the Tree sub-module...");
tem   = FT(25);
wind  = FT(2);
width = FT(0.05);
@btime radiative_conductance($tem);
@btime black_body_emittance($tem);
@btime boundary_layer_conductance($wind, $width);
=#




# Benchmarking functions
tree   = Plants.Yujie111{FT}();
envir  = AirLayer{FT}();
c3_set = C3CLM(FT);
zenith = FT(30);
r_all  = FT(1000);
Plants.Yujie111EvaluateModel(tree, c3_set, envir, zenith, r_all);

