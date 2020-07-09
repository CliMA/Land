# include("./scripts/safety_efficiency_tradeoff/tradeoff_3_generate_optima_scatter_thread.jl")




# function to run all required result all at once
function RunTradeoffOptimaScatter(weather, years, filename, d_lati, d_long, d_alti, selections)
    # initialize the node
    @everywhere max_vmax = 100.0
    @everywhere opt_node = Yujie111Init()
    @everywhere opt_node.k_root *= 0.2
    @everywhere opt_node.k_stem *= 0.2
    @everywhere opt_node.k_sla  *= 0.2
    @everywhere opt_node.d_lati  = $d_lati
    @everywhere opt_node.d_long  = $d_long
    @everywhere opt_node.d_alti  = $d_alti
    @everywhere Yujie111UpdateSoilFromSWC(opt_node, 1.0)
    @everywhere Yujie111UpdateLeaf(opt_node, 1000.0, 70.0)

    # define the thread function
    @everywhere weat_years = CSV.read($weather)

    @everywhere function GetOptimalBCKLVForYear(param)
        nth       = param[1]
        year      = param[2]
        selection = param[3]
        new_b     = param[4]
        disp      = param[5]
        tmp_node  = deepcopy(opt_node)
        tmp_node.b_root = new_b
        tmp_node.b_stem = new_b
        tmp_node.b_leaf = new_b
        ini_laba  = opt_node.laba
        ini_vmax  = opt_node.v_max
        println("Start the ", nth, " th simulation!")
        println("Start year ", year, "!")
        weat_mask = (weat_years.Year .== year)
        weat_year = weat_years[weat_mask,:]
        result = Yujie111GetOptimalBCKLV(tmp_node, weat_year, selection, max_vmax, disp)
        return result
    end

    # run the simulations
    params = []
    count  = 0
    newbs  = 1.0:0.5:10.0
    total  = length(newbs) * length(years)
    for newb in newbs
        for year in years
            count += 1
            if count<=5 || count>=total-5
                push!(params, [count year selections newb true ])
            else
                push!(params, [count year selections newb false])
            end
        end
    end
    result = pmap(GetOptimalBCKLVForYear, params)
    writedlm(filename, result)
end




# main function
list_loca = CSV.read( homedir() * "/Data/USAForest20Sites/location.csv" )
list_site = list_loca.Site
list_lati = list_loca.Latitude
list_long = list_loca.Longitude
list_alti = list_loca.Elevation

for i in [15]#[4,7,8,15,20]
    println(i)
    
    site = list_site[i]
    lati = list_lati[i]
    long = list_long[i]
    alti = list_alti[i]

    weat_his_amb  =  homedir() * "/Data/USAForest20Sites/" * site * "/" * site * "_histo_sim0_GS.csv"

    RunTradeoffOptimaScatter(weat_his_amb, 1971:2005, "./data/" * site * "/opt_bk_scatter.txt", lati, long, alti, "KLV")
end