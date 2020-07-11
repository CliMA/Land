# include("src/test.jl")

using BenchmarkTools
using CSV
using DataFrames
using Photosynthesis
using Plants

PL = Plants
FT = Float64

node  = PL.SPACSimple{FT}();
photo = C3CLM(FT);
zeni  = FT(30);
rall  = FT(1000);

#weat_hist = DataFrame!(CSV.File("./data/Flagstaff_histo_sim0_GS.csv"))
weat_hist = DataFrame!(CSV.File("/home/wyujie/Data/USAForest20Sites/Flagstaff/Flagstaff_histo_sim0_GS.csv"))
weat_mask = (weat_hist.Year .== 2005);
weat_2005 = weat_hist[weat_mask,:];


PL.Yujie111EvaluateModel(node, photo, zeni, rall);
PL.Yujie111GainRiskMatrix(node, photo, zeni, rall);
PL.Yujie111GetAnnualCiCa(node, photo, weat_2005);
PL.Yujie111GetAnnualProfit(node, photo, weat_2005, FT(100));
PL.Yujie111GetOptimaldAdE(node, photo, zeni, rall);
PL.Yujie111GetOptimalFs(node, photo, zeni, rall);
#PL.Yujie111GetOptimalFsMap(node, photo, zeni, rall);
PL.Yujie111GetOptimalInvestment(node, photo, weat_2005);

#@btime PL.Yujie111EvaluateModel(node, photo, zeni, rall);





#=
include("./include_leaf_rubisco.jl")

using DelimitedFiles


palisade = YujiePalisade1DInit()
list_order = Palisade1DAllocateRubiscoNS(palisade, 200, true)
writedlm("./allocation_dynamics.txt", list_order)

palisade = YujiePalisade1DInit()
list_unit = 0
list_anet = 0
last_anet = -9999.9
for i in 1:200
    println(i)
    global last_anet,list_unit,list_anet
    Palisade1DAllocateRubiscoNS(palisade, 1)
    anet = Palisade1DGetAssimilation(palisade)
    #=
    if anet>last_anet
        list_unit = [list_unit; i   ]
        list_anet = [list_anet; anet]
        last_anet = anet
    else
        break
    end
    =#
    list_unit = [list_unit; i   ]
    list_anet = [list_anet; anet]
end
writedlm("./allocation_dynamics_anet.txt", [list_unit list_anet])
=#

#= test the chlorolast re-distribution
palisade = YujiePalisade1DInit()
Palisade1DAllocateRubisco(palisade, 100)

a = Palisade1DGetAssimilation(palisade)
println(a)

palisade.co2 = 100
a = Palisade1DGetAssimilation(palisade)
println(a)

list_order = []
Palisade1DReallocateRubisco(palisade, list_order)
a = Palisade1DGetAssimilation(palisade)
println(a)
=#

#=
mm = Palisade1DOptimizeleaf()
maxl = 0
for i in 1:length(mm)
    global maxl
    if length(mm[i])>maxl
        maxl = length(mm[i])
    end
end

nn = zeros(length(mm), maxl)
for i in 1:length(mm)
    global maxl
    for j in 1:length(mm[i])
        nn[i,j] = max(mm[i][j], 0)
    end
end

lai = 0.1:0.1:length(mm)*0.1
rub = 1:maxl

println(nn)
=#

#palisade = YujiePalisade1DInit()
#Palisade1DAllocateRubisco(palisade, 120, true)
#Palisade1DUpdateACc(palisade)

#mat = [palisade.r_list palisade.l_list palisade.c_list palisade.a_list]
#writedlm("/home/yujie/allocation_example.txt", mat)

#list_order = Palisade1DAllocateRubisco(palisade, 200, true)
#writedlm("/home/yujie/allocation_dynamics.txt", list_order)

#Palisade1DMaximumVmaxTemperature(palisade)
#Palisade1DMaximumVmaxCO2(palisade)
#Palisade1DMaximumVmaxPAR(palisade)
#Palisade1DMaximumVmaxGwlmm(palisade)
#Palisade1DMaximumVmaxUnit(palisade)

#=
list_rm = []
list_rr = []
list_am = []

cell = YujiePalisade1DInit()

max_rub = 0
for rub in 5:5:300
    global max_rub
    cell.par = 1600
    cell.co2 = 30.0
    list_c = []
    list_a = []
    list_t = []
    Palisade1DAllocateRubisco(cell, 5)
    a = Palisade1DGetAssimilation(cell)
    rur = sum(cell.r_list)
    if rur>max_rub
        max_rub = rur
    else
        break
    end
    push!(list_rm, rub)
    push!(list_rr, rur)
    push!(list_am, a  )
    file_name = "/home/yujie/Documents/rub-" * string(Int(rur)) * "-aci.txt"
    println('"', file_name, '"')
    for co2 in 5:5:200
        cell.co2 = co2
        a = Palisade1DGetAssimilation(cell)
        push!(list_c, co2     )
        push!(list_a, a       )
        push!(list_t, cell.tem)
    end
    writedlm(file_name, [list_c list_a list_t])
end
=#









#=
function TestNode111()
    node = YujieNode111Init()
    YujieNode111EvaluateNightE(node)
    #=
    node = YujieNode111Init()
    node.b_root = 1.879
    node.c_root = 2.396
    node.b_stem = 2.238
    node.c_stem = 9.380
    node.b_leaf = 1.897
    node.c_leaf = 2.203
    node.k_root = 2228.0
    node.k_stem = 4926.7
    node.k_leaf = 5424.3
    YujieNode111UpdateLeaf(node, 4758.5, 61.74)
    YujieNode111UpdateSoilFromP(node, 0.8)
    YujieNode111UpdateLegacy(node, 0.0)
    tmp_mat = [node.l_root; node.l_stem; node.l_leaf]
    for vpd in [1.0, 1.5, 2.0, 2.5, 3.0]
    #for co2 in [20, 30, 40, 60, 80]
    #for ps in [0.0, 0.4, 0.8, 1.2, 1.6]
        tmp_node = deepcopy(node)
        #YujieNode111UpdateSoilFromP(tmp_node, ps)
        tmp_result = YujieNode111GetOptima(tmp_node, 25.0, 750.0, vpd, 40.0, 1000.0)
        #tmp_result = YujieNode111GetOptima(tmp_node, 25.0, 750.0, 1.5, co2, 1000.0)
        #tmp_result = YujieNode111GetOptima(tmp_node, 25.0, 750.0, 1.5, 40.0, 1000.0)
        println(tmp_result)
        flow = tmp_result["f"]
        YujieNode111UpdateLegacy(tmp_node, flow)
        tmp_list = [tmp_node.l_root; tmp_node.l_stem; tmp_node.l_leaf]
        tmp_mat  = [tmp_mat tmp_list]
    end
    writedlm("/home/yujie/example_legacy.txt", tmp_mat)
    =#
end

function TestNode112()
    node = YujieNode112Init(0.5, 0.5)
    node.p_soil = 0.5
    YujieNode112PlotGainRisk(node)
end
=#