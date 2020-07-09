# include("./scripts/safety_efficiency_tradeoff/tradeoff_4_make_sensitivity_analysis_thread.jl")




# function to run all required result all at once
function RunTradeoffSensitivityAnalysis(weather, years, filename, d_lati, d_long, d_alti, factor_to_vary, ratio_to_vary, selection)
    # initialize the node
    @everywhere max_vmax = 100.0
    @everywhere opt_node = Yujie111Init()
    @everywhere opt_node.k_root *= 0.2
    @everywhere opt_node.k_stem *= 0.2
    @everywhere opt_node.k_sla  *= 0.2
    @everywhere opt_node.d_lati = $d_lati
    @everywhere opt_node.d_long = $d_long
    @everywhere opt_node.d_alti = $d_alti
    @everywhere Yujie111UpdateSoilFromSWC(opt_node, 1.0)
    @everywhere Yujie111UpdateLeaf(opt_node, 1000.0, 70.0)

    # define the thread function
    @everywhere weat_years = CSV.read($weather)
    if factor_to_vary=="cv"
        @everywhere opt_node.c_vmax *= $ratio_to_vary
    elseif factor_to_vary=="cc"
        @everywhere opt_node.c_cons *= $ratio_to_vary
    elseif factor_to_vary=="ge"
        @everywhere opt_node.g_mul  *= $ratio_to_vary
    elseif factor_to_vary=="gm"
        @everywhere opt_node.g_sla  *= $ratio_to_vary
    elseif factor_to_vary=="ga"
        @everywhere opt_node.gaba   *= $ratio_to_vary
    elseif factor_to_vary=="sd"
        @everywhere opt_node.h_root *= $ratio_to_vary
        @everywhere opt_node.h_soil *= $ratio_to_vary
    elseif factor_to_vary=="ca"
        @everywhere  weat_years.CO2 *= $ratio_to_vary
    elseif factor_to_vary=="rh"
        @everywhere weat_years.D    *= $ratio_to_vary
    elseif factor_to_vary=="ta"
        @everywhere temp_rh_multip1  = GetSaturatedVaporPressureList(weat_years.Tair)
        @everywhere weat_years.Tair  = weat_years.Tair .+ $ratio_to_vary
        @everywhere temp_rh_multip2  = GetSaturatedVaporPressureList(weat_years.Tair)
        @everywhere weat_years.D     = weat_years.D .* temp_rh_multip2 ./ temp_rh_multip1
    end

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
                push!(params, [count year selection newb true ])
            else
                push!(params, [count year selection newb false])
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

for i in [7]#[4,7,8,15,20]
    println(i)
    
    site = list_site[i]
    lati = list_lati[i]
    long = list_long[i]
    alti = list_alti[i]

    weat_his_amb  =  homedir() * "/Data/USAForest20Sites/" * site * "/" * site * "_histo_sim0_GS.csv"

    # run default settings first
    RunTradeoffSensitivityAnalysis(weat_his_amb, 1971:2005, "./data/" * site * "/optk_def.txt", lati, long, alti, "def", 1.0, "KLV")

    # sensitivity analysis
    for lab in ["ca", "cc", "cv", "ga", "gm", "rh", "sd", "ta"]
        println("Nowing run sensitivity to ", lab)
        if lab=="rh"
            for ratio_to_vary in [0.8,1.2]
                println("Now running its ", ratio_to_vary, "X the default")
                file_optk = "./data/" * site * "/optk_" * lab * "_" * string(ratio_to_vary) * ".txt"
                RunTradeoffSensitivityAnalysis(weat_his_amb, 1971:2005, file_optk, lati, long, alti, lab, ratio_to_vary, "KLV")
            end
        elseif lab=="ta"
            for ratio_to_vary in [-2,2]
                if ratio_to_vary<0
                    println("Now running its default", ratio_to_vary)
                    file_optk = "./data/" * site * "/optk_" * lab * "_" * string(ratio_to_vary) * ".txt"
                else
                    println("Now running its default+", ratio_to_vary)
                    file_optk = "./data/" * site * "/optk_" * lab * "_+" * string(ratio_to_vary) * ".txt"
                end
                RunTradeoffSensitivityAnalysis(weat_his_amb, 1971:2005, file_optk, lati, long, alti, lab, ratio_to_vary, "KLV")
            end
        else
            for ratio_to_vary in [0.5,1.5]
                println("Now running its ", ratio_to_vary, "X the default")
                file_optk = "./data/" * site * "/optk_" * lab * "_" * string(ratio_to_vary) * ".txt"
                RunTradeoffSensitivityAnalysis(weat_his_amb, 1971:2005, file_optk, lati, long, alti, lab, ratio_to_vary, "KLV")
            end
        end
    end
end