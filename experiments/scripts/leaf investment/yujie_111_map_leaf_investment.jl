# include("./scripts/leaf investment/yujie_111_map_leaf_investment.jl")

# make node for year 2005
node = Yujie111Init()
node.d_lati =  35.198284
node.d_long = -111.651299
node.d_alti =  2105.0
Yujie111UpdateLeaf(node, 1193.0, 88.9)
Yujie111UpdateSoilFromSWC(node, 1.0)

#weat_hist = DataFrame!(CSV.File(homedir() * "/Data/USAForest20Sites/Flagstaff/Flagstaff_histo_sim0_GS.csv"))
weat_hist = DataFrame!(CSV.File("./data/Flagstaff_histo_sim0_GS.csv"))
weat_mask = (weat_hist.Year .== 2005)
weat_2005 = weat_hist[weat_mask,:]

mapp = []
for laba in 200:200:5000
    global mapp
    list_p = []
    for vmax in 10:5:100
        tmp_node = deepcopy(node)
        Yujie111UpdateLeaf(tmp_node, laba, vmax)
        prof = Yujie111GetAnnualProfit(tmp_node, weat_2005, 100.0)
        println("LaBa: ", laba, "\tVmax: ", vmax, "\tProfit: ", prof)
        if list_p==[]
            list_p = prof
        else
            list_p = [list_p prof]
        end
    end
    println(list_p)
    if mapp==[]
        mapp = list_p
    else
        mapp =[mapp; list_p]
    end
end

writedlm("./optimization_map_2005.txt", mapp)





#=
# make node for year 2097
node = Yujie111Init()
node.d_lati =  35.198284
node.d_long = -111.651299
node.d_alti =  2105.0
Yujie111UpdateLeaf(node, 2229.0, 69.25)
Yujie111UpdateSoilFromSWC(node, 1.0)

weat_futu = CSV.read("./data/Flagstaff/rcp85_30_ele.csv")
weat_mask = (weat_futu.Year .== 2097)
weat_2097 = weat_futu[weat_mask,:]

mapp = []
for laba in 200:200:5000
    global mapp
    list_p = []
    for vmax in 10:5:100
        tmp_node = deepcopy(node)
        Yujie111UpdateLeaf(tmp_node, laba, vmax)
        prof = Yujie111GetAnnualProfit(tmp_node, weat_2097, 100.0)
        println("LaBa: ", laba, "\tVmax: ", vmax, "\tProfit: ", prof)
        if list_p==[]
            list_p = prof
        else
            list_p = [list_p prof]
        end
    end
    println(list_p)
    if mapp==[]
        mapp = list_p
    else
        mapp =[mapp; list_p]
    end
end

writedlm("./optimization_map_2097.txt", mapp)
=#