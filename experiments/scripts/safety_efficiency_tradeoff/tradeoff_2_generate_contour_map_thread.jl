# include("./scripts/safety_efficiency_tradeoff/tradeoff_2_generate_contour_map_thread.jl")




# define the map function to obtain matrix with threaded function
function GenerateBCKContourMap(weather, year, filename, d_lati, d_long, d_alti)
    # Define a default node and default weather
    @everywhere node_default = Yujie111Init()
    @everywhere node_default.d_lati = $d_lati
    @everywhere node_default.d_long = $d_long
    @everywhere node_default.d_alti = $d_alti
    @everywhere Yujie111UpdateSoilFromSWC(node_default, 1.0)
    @everywhere weat_years = CSV.read($weather)

    # define threading function
    @everywhere function GenerateGSCPFromParameter(param)
        year = $year
        coun = param[1]
        krat = param[2]
        bval = param[3]
        println("\t", year, "\t", coun, "\t", krat, "\t", bval)
        weat_mask = (weat_years.Year .== year)
        weat_year = weat_years[weat_mask,:]

        node_opt = deepcopy(node_default)
        node_opt.k_root *= krat
        node_opt.k_stem *= krat
        node_opt.k_sla  *= krat
        node_opt.b_root  = bval
        node_opt.b_stem  = bval
        node_opt.b_leaf  = bval

        optl,optv,optp = Yujie111GetOptimalInvestment(node_opt, weat_year, 2000.0, 60.0, 100.0, false)

        return [krat bval optl optv optp]
    end

    # generate parameter list
    list_param = []
    count      = 0
    for krat_tmp in 0.1:0.1:2.0
        for bval_tmp in 0.5:0.25:10.0
            count += 1
            push!(list_param, [count krat_tmp bval_tmp])
        end
    end
    result = pmap(GenerateGSCPFromParameter, list_param)

    # output the results
    CSV.write(filename, DataFrame( vcat(result...) ), header=["KRAT","BVAL","LABA","VMAX","GSCP"])
end




# function to run through several sites
list_loca = CSV.read( homedir() * "/Data/USAForest20Sites/location.csv" )
list_site = list_loca.Site
list_lati = list_loca.Latitude
list_long = list_loca.Longitude
list_alti = list_loca.Elevation

# Durango      4
# Flagstaff    7 done
# Hattiesburg  8
# Northway    15
# Trinity     20

for i in [7]#[4,7,8,15,20]
    println(i)
    
    site = list_site[i]
    lati = list_lati[i]
    long = list_long[i]
    alti = list_alti[i]

    weat_his_amb =  homedir() * "/Data/USAForest20Sites/" * site * "/" * site * "_histo_sim0_GS.csv"
    for year in 1971:2004#2005
        file_to_write = "./data/" * site * "_contour_" * string(year) * ".csv"
        GenerateBCKContourMap(weat_his_amb, year, file_to_write, lati, long, alti)
    end
end