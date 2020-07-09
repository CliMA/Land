# include("./scripts/leaf investment/yujie_111_optima_sensitivity.jl")

# for current ambient
function GenerateOptimaPTh(weather, years, filename, d_lati, d_long, d_alti, label)
    # initialize the node
    max_vmax = 100.0
    ini_vmax = 70.0
    ini_laba = 1500.0
    opt_node = Yujie111Init()
    opt_node.d_lati = d_lati
    opt_node.d_long = d_long
    opt_node.d_alti = d_alti
    if label=="B"
        opt_node.b_root *= 2.0
        opt_node.b_stem *= 2.0
        opt_node.b_leaf *= 2.0
    elseif label=="C"
        opt_node.c_root *= 2.0
        opt_node.c_stem *= 2.0
        opt_node.c_leaf *= 2.0
    elseif label=="K"
        opt_node.k_root *= 2.0
        opt_node.k_stem *= 2.0
        opt_node.k_sla  *= 2.0
    elseif label=="KW"
        opt_node.k_root *= 2.0
        opt_node.k_stem *= 2.0
    elseif label=="CV"
        opt_node.c_vmax *= 2.0
    elseif label=="CC"
        opt_node.c_cons *= 2.0
    elseif label=="GCE"
        opt_node.g_mul  *= 2.0
    elseif label=="GABA"
        opt_node.gaba   *= 2.0
    end
    Yujie111UpdateSoilFromSWC(opt_node, 1.0)
    Yujie111UpdateLeaf(opt_node, 2000.0, 50.0)

    # define the thread function
    weat_years = DataFrame!(CSV.File(weather))
    weat_mask = (weat_years.Year .== 2005)
    weat_year = weat_years[weat_mask,:]
    opt_laba,opt_vmax,opt_prof = Yujie111GetOptimalInvestment(opt_node, weat_year, ini_laba, ini_vmax, max_vmax)
    Yujie111UpdateLeaf(opt_node, opt_laba, opt_vmax)
    opt_cica = Yujie111GetAnnualCiCa(opt_node, weat_year)
    println("Finished label ", label, "!\n\n")
    return [label opt_laba opt_vmax opt_prof opt_cica[1] opt_cica[2]]
end

function GenerateOptimaNTh(weather, years, filename, d_lati, d_long, d_alti, label)
    # initialize the node
    max_vmax = 100.0
    ini_vmax = 70.0
    ini_laba = 1500.0
    opt_node = Yujie111Init()
    opt_node.d_lati = d_lati
    opt_node.d_long = d_long
    opt_node.d_alti = d_alti
    if label=="B"
        opt_node.b_root *= 0.5
        opt_node.b_stem *= 0.5
        opt_node.b_leaf *= 0.5
    elseif label=="C"
        opt_node.c_root *= 0.5
        opt_node.c_stem *= 0.5
        opt_node.c_leaf *= 0.5
    elseif label=="K"
        opt_node.k_root *= 0.5
        opt_node.k_stem *= 0.5
        opt_node.k_sla  *= 0.5
    elseif label=="KW"
        opt_node.k_root *= 0.5
        opt_node.k_stem *= 0.5
    elseif label=="CV"
        opt_node.c_vmax *= 0.5
    elseif label=="CC"
        opt_node.c_cons *= 0.5
    elseif label=="GCE"
        opt_node.g_mul  *= 0.5
    elseif label=="GABA"
        opt_node.gaba   *= 0.5
    end
    Yujie111UpdateSoilFromSWC(opt_node, 1.0)
    Yujie111UpdateLeaf(opt_node, 2000.0, 50.0)

    # define the thread function
    weat_years = DataFrame!(CSV.File(weather))
    weat_mask = (weat_years.Year .== 2005)
    weat_year = weat_years[weat_mask,:]
    opt_laba,opt_vmax,opt_prof = Yujie111GetOptimalInvestment(opt_node, weat_year, ini_laba, ini_vmax, max_vmax)
    Yujie111UpdateLeaf(opt_node, opt_laba, opt_vmax)
    opt_cica = Yujie111GetAnnualCiCa(opt_node, weat_year)
    println("Finished label ", label, "!\n\n")
    return [label opt_laba opt_vmax opt_prof opt_cica[1] opt_cica[2]]
end




# function to run through several sites
list_site = ["Flagstaff", "Hattiesburg", "Missoula" , "Trinity" ]
list_lati = [35.198284  , 31.32712     , 46.965260  , 30.941690 ]
list_long = [-111.651299, -89.290337   , -109.533691, -95.377121]
list_alti = [2105.0     , 50.0         , 1344.0     , 68.0      ]

for i in 1:1
    site = list_site[i]
    lati = list_lati[i]
    long = list_long[i]
    alti = list_alti[i]

    weat_his_amb = "./data/" * site * "/histo_30_amb.csv"

    year_his = 1971:2005
    
    list_result = []
    for label in ["DEF" "B" "C" "K" "CV" "CC" "GCE" "GABA"]
        tmp_re = GenerateOptimaNTh(weat_his_amb, year_his, "./data/" * site * "/opt_his_amb.txt", lati, long, alti, label)
        if list_result==[]
            list_result = tmp_re
        else
            list_result = [list_result; tmp_re]
        end
    end

    println(list_result)
    writedlm("./data/" * site * "/sensitivity_analysis_n.txt", list_result)
    
    list_result = []
    for label in ["DEF" "B" "C" "K" "CV" "CC" "GCE" "GABA"]
        tmp_re = GenerateOptimaPTh(weat_his_amb, year_his, "./data/" * site * "/opt_his_amb.txt", lati, long, alti, label)
        if list_result==[]
            list_result = tmp_re
        else
            list_result = [list_result; tmp_re]
        end
    end

    println(list_result)
    writedlm("./data/" * site * "/sensitivity_analysis_p.txt", list_result)
end