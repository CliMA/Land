# include("scripts/yujie_111_price_of_water_2.jl")
file_list = ["./price_of_water_1.txt",
             "./price_of_water_2.txt",
             "./price_of_water_3.txt",
             "./price_of_water_4.txt",
             "./price_of_water_5.txt"]


# ca response
ca_list = [20.0, 30.0, 40.0, 60.0, 80.0]
for i in 1:5
    node   = Yujie111Init();
    Yujie111UpdateLeaf(node, 750.0, 100.0)
    Yujie111UpdateSoilFromP(node, 0.0);
    envir  = [25.0 ca_list[i] 1.5 2.0];
    zenith = 30.0
    r_all  = 1000.0
    matrix = Yujie111EvaluateModel(node, envir, zenith, r_all)
    writedlm(file_list[i], matrix)
end


#=
# d response
d_list = [0.5, 1.0, 1.5, 2.0, 2.5]
for i in 1:5
    node   = Yujie111Init();
    Yujie111UpdateLeaf(node, 750.0, 100.0)
    Yujie111UpdateSoilFromP(node, 0.0);
    envir  = [25.0 40.0 d_list[i] 2.0];
    zenith = 30.0
    r_all  = 1000.0
    matrix = Yujie111EvaluateModel(node, envir, zenith, r_all)
    writedlm(file_list[i], matrix)
end
=#

#=
# ps response
p_list = [0.0, 0.5, 1.0, 1.5, 2.0]
for i in 1:5
    node   = Yujie111Init();
    Yujie111UpdateLeaf(node, 750.0, 100.0)
    Yujie111UpdateSoilFromP(node, p_list[i]);
    envir  = [25.0 40.0 1.5 2.0];
    zenith = 30.0
    r_all  = 1000.0
    matrix = Yujie111EvaluateModel(node, envir, zenith, r_all)
    writedlm(file_list[i], matrix)
end
=#
