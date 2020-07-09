# include("./scripts/leaf investment/yujie_111_map_leaf_cica.jl")

# make node for year 2005
node = Yujie111Init()
node.d_lati =  35.198284
node.d_long = -111.651299
node.d_alti =  2105.0
Yujie111UpdateLeaf(node, 1446.0, 84.1)
Yujie111UpdateSoilFromSWC(node, 1.0)

#weat_hist = DataFrame!(CSV.File(homedir() * "/Data/USAForest20Sites/Flagstaff/Flagstaff_histo_sim0_GS.csv"))
weat_hist = DataFrame!(CSV.File("./data/Flagstaff_histo_sim0_GS.csv"))
weat_2005 = weat_hist[weat_mask,:]

mapp = []
for laba in 200:200:5000
    global mapp
    list_p = []
    for vmax in 10:5:100
        tmp_node = deepcopy(node)
        Yujie111UpdateLeaf(tmp_node, laba, vmax)
        cica = Yujie111GetAnnualCiCa(tmp_node, weat_2005)
        println("LaBa: ", laba, "\tVmax: ", vmax, "\tCi/Ca: ", cica[1])
        if list_p==[]
            list_p = cica[1]
        else
            list_p = [list_p cica[1]]
        end
    end
    println(list_p)
    if mapp==[]
        mapp = list_p
    else
        mapp =[mapp; list_p]
    end
end

writedlm("./cica_map_2005.txt", mapp)





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
        cica = Yujie111GetAnnualCiCa(tmp_node, weat_2097)
        println("LaBa: ", laba, "\tVmax: ", vmax, "\tCi/Ca: ", cica)
        if list_p==[]
            list_p = cica
        else
            list_p = [list_p cica]
        end
    end
    println(list_p)
    if mapp==[]
        mapp = list_p
    else
        mapp =[mapp; list_p]
    end
end

writedlm("./cica_map_2097.txt", mapp)
=#
