# include("./scripts/leaf investment/yujie_111_map_optima_shift.jl")

# for current ambient
function GenerateOptima(weather, years, filename, d_lati, d_long, d_alti)
    # initialize the node
    @everywhere max_vmax = 100.0
    @everywhere opt_node = Yujie111Init()
    @everywhere opt_node.d_lati = $d_lati
    @everywhere opt_node.d_long = $d_long
    @everywhere opt_node.d_alti = $d_alti
    @everywhere Yujie111UpdateSoilFromSWC(opt_node, 1.0)
    @everywhere Yujie111UpdateLeaf(opt_node, 1000.0, 70.0)

    # define the thread function
    @everywhere weat_years = DataFrame!(CSV.File($weather))
    @everywhere function GetOptimaForYear(year)
        tmp_node = deepcopy(opt_node)
        ini_laba = opt_node.laba
        ini_vmax = opt_node.v_max
        println("Start year ", year, "!")
        weat_mask = (weat_years.Year .== year)
        weat_year = weat_years[weat_mask,:]
        opt_laba,opt_vmax,opt_prof = Yujie111GetOptimalInvestment(tmp_node, weat_year, ini_laba, ini_vmax, max_vmax)
        Yujie111UpdateLeaf(opt_node, opt_laba, opt_vmax)
        println("Finished year ", year, "!\n\n")
        return [opt_laba opt_vmax opt_prof]
    end

    # run the simulations
    list_year = years
    result = pmap(GetOptimaForYear, list_year)
    writedlm(filename, result)
end




# function to run through several sites
list_loca = DataFrame!(CSV.File("./data/location.csv"))
list_site = list_loca.Site
list_lati = list_loca.Latitude
list_long = list_loca.Longitude
list_alti = list_loca.Elevation

# Flagstaff 7
# Hattiesburg 8
# Missoula 13
# Trinity 20
# Finished 1 7 8 13 20

for i in [2 3 4 5 6 9 10 11 12 14 15 16 17 18 19]
    print(i)
    
    site = list_site[i]
    lati = list_lati[i]
    long = list_long[i]
    alti = list_alti[i]

    weat_his_amb = "./data/" * site * "/histo_30_amb.csv"
    weat_his_ele = "./data/" * site * "/histo_30_ele.csv"
    weat_fut_amb = "./data/" * site * "/rcp85_30_amb.csv"
    weat_fut_ele = "./data/" * site * "/rcp85_30_ele.csv"

    year_his = 1971:2005
    year_fut = 2065:2099

    GenerateOptima(weat_his_amb, year_his, "./data/" * site * "/opt_his_amb.txt", lati, long, alti)
    GenerateOptima(weat_his_ele, year_his, "./data/" * site * "/opt_his_ele.txt", lati, long, alti)
    GenerateOptima(weat_fut_amb, year_fut, "./data/" * site * "/opt_fut_amb.txt", lati, long, alti)
    GenerateOptima(weat_fut_ele, year_fut, "./data/" * site * "/opt_fut_ele.txt", lati, long, alti)
end