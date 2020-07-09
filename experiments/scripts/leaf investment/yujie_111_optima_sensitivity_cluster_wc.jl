# include("./scripts/leaf investment/yujie_111_optima_sensitivity_cluster_wc.jl")

# for current ambient
function GenerateOptima(weather, years, filename, d_lati, d_long, d_alti, factor_to_vary, ratio_to_vary)
    # initialize the node
    @everywhere max_vmax = 100.0
    @everywhere opt_node = Yujie111Init()
    @everywhere opt_node.c_root = 3.0
    @everywhere opt_node.c_stem = 3.0
    @everywhere opt_node.c_leaf = 3.0

    @everywhere opt_node.d_lati = $d_lati
    @everywhere opt_node.d_long = $d_long
    @everywhere opt_node.d_alti = $d_alti
    @everywhere Yujie111UpdateSoilFromSWC(opt_node, 1.0)
    @everywhere Yujie111UpdateLeaf(opt_node, 1000.0, 70.0)

    # define the thread function
    @everywhere weat_years = DataFrame!(CSV.File($weather))

    if factor_to_vary=="wb"
        @everywhere opt_node.b_root *= $ratio_to_vary
        @everywhere opt_node.b_stem *= $ratio_to_vary
        @everywhere opt_node.b_leaf *= $ratio_to_vary
    elseif factor_to_vary=="wc"
        @everywhere opt_node.c_root *= $ratio_to_vary
        @everywhere opt_node.c_stem *= $ratio_to_vary
        @everywhere opt_node.c_leaf *= $ratio_to_vary
    elseif factor_to_vary=="wk"
        @everywhere opt_node.k_root *= $ratio_to_vary
        @everywhere opt_node.k_stem *= $ratio_to_vary
        @everywhere opt_node.k_sla  *= $ratio_to_vary
    elseif factor_to_vary=="kw"
        @everywhere opt_node.k_root *= $ratio_to_vary
        @everywhere opt_node.k_stem *= $ratio_to_vary
    elseif factor_to_vary=="kl"
        @everywhere opt_node.k_sla  *= $ratio_to_vary
    elseif factor_to_vary=="cv"
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

    @everywhere function GetOptimaForYear(year)
        tmp_node = deepcopy(opt_node)
        ini_laba = opt_node.laba
        ini_vmax = opt_node.v_max
        println("Start year ", year, "!")
        weat_mask = (weat_years.Year .== year)
        weat_year = weat_years[weat_mask,:]
        opt_laba,opt_vmax,opt_prof = Yujie111GetOptimalInvestment(tmp_node, weat_year, ini_laba, ini_vmax, max_vmax)
        Yujie111UpdateLeaf(opt_node, opt_laba, opt_vmax)
        opt_cica = Yujie111GetAnnualCiCa(opt_node, weat_year)
        println("Finished year ", year, "!\n\n")
        return [opt_laba opt_vmax opt_prof opt_cica[1]]
    end

    # run the simulations
    list_year = years
    result = pmap(GetOptimaForYear, list_year)
    writedlm(filename, result)
end




# function to run through several sites
#loca_file = homedir() * "\\Box Sync\\Research_Ongoing\\SPAC - LeafOptimization\\shifts\\location.csv"
loca_file = homedir() * "/Data/USAForest20Sites/location.csv"
list_loca = DataFrame!(CSV.File(loca_file))
list_site = list_loca.Site
list_lati = list_loca.Latitude
list_long = list_loca.Longitude
list_alti = list_loca.Elevation

# Durango      4
# Flagstaff    7
# Hattiesburg  8
# Trinity     20

for i in [7]
    println(i)
    
    site = list_site[i]
    lati = list_lati[i]
    long = list_long[i]
    alti = list_alti[i]

    #weat_his_amb =  homedir() * "\\Box Sync\\Research_Ongoing\\SPAC - LeafOptimization\\shifts\\" * site * "\\" * site * "_histo_sim0_GS.csv"
    weat_his_amb =  homedir() * "/Data/USAForest20Sites/" * site * "/" * site * "_histo_sim0_GS.csv"

    year_his = 1971:2005

    GenerateOptima(weat_his_amb, year_his, "./data/" * site * "/c3/amb_def.txt", lati, long, alti, "wb", 1.0)
    #for lab in ["ca", "cc", "cv", "ga", "ge", "gm", "kl", "kw", "rh", "sd", "ta", "wb", "wc", "wk"]
    for lab in ["wb"]#["ca", "cc", "cv", "ga", "gm", "kl", "kw", "rh", "sd", "ta", "wb", "wc", "wk"]
        println("Nowing run sensitivity to ", lab)
        if lab=="ca"
            for ratio_to_vary in [0.5,0.75,1.5,2.0]
                println("Now running its ", ratio_to_vary, "X the default")
                file_to_write = "./data/" * site * "/c3/amb_" * lab * "_" * string(ratio_to_vary) * ".txt"
                GenerateOptima(weat_his_amb, year_his, file_to_write, lati, long, alti, lab, ratio_to_vary)
            end
        elseif lab=="rh"
            for ratio_to_vary in [0.6,0.8,1.2,1.4]
                println("Now running its ", ratio_to_vary, "X the default")
                file_to_write = "./data/" * site * "/c3/amb_" * lab * "_" * string(ratio_to_vary) * ".txt"
                GenerateOptima(weat_his_amb, year_his, file_to_write, lati, long, alti, lab, ratio_to_vary)
            end
        elseif lab=="ta"
            for ratio_to_vary in [-2,2,4,6]
                if ratio_to_vary<0
                    println("Now running its default", ratio_to_vary)
                    file_to_write = "./data/" * site * "/c3/amb_" * lab * "_" * string(ratio_to_vary) * ".txt"
                else
                    println("Now running its default+", ratio_to_vary)
                    file_to_write = "./data/" * site * "/c3/amb_" * lab * "_+" * string(ratio_to_vary) * ".txt"
                end
                GenerateOptima(weat_his_amb, year_his, file_to_write, lati, long, alti, lab, ratio_to_vary)
            end
        else
            for ratio_to_vary in [0.2,0.6,1.5,2.0]
                println("Now running its ", ratio_to_vary, "X the default")
                file_to_write = "./data/" * site * "/c3/amb_" * lab * "_" * string(ratio_to_vary) * ".txt"
                GenerateOptima(weat_his_amb, year_his, file_to_write, lati, long, alti, lab, ratio_to_vary)
            end
        end
    end
end