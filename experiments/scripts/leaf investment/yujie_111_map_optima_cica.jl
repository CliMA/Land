# include("./scripts/leaf investment/yujie_111_map_optima_cica.jl")

# function to get cica
function GetCiCas(weather, years, optima, d_lati, d_long, d_alti)
    tmp_node = Yujie111Init()
    tmp_node.d_lati = d_lati
    tmp_node.d_long = d_long
    tmp_node.d_alti = d_alti
    list_cica = []
    for i in 1:length(years)
        Yujie111InitializeLegacy(tmp_node)
        Yujie111UpdateSoilFromSWC(tmp_node, 1.0)
        Yujie111UpdateLeaf(tmp_node, optima[i,1], optima[i,2])
        weat_mask  = (weather.Year .== years[i])
        weat_year  = weather[weat_mask,:]
        tmp_result = Yujie111GetAnnualCiCa(tmp_node, weat_year)
        println(years[i], "\t", tmp_result)
        if list_cica==[]
            list_cica = tmp_result
        else
            list_cica = [list_cica; tmp_result]
        end
    end
    return list_cica
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
# Finished 7 8 13 20

for i in [2 3 4 5 6 9 10 11 12 14 15 16 17 18 19]
    println(i)
    
    site = list_site[i]
    lati = list_lati[i]
    long = list_long[i]
    alti = list_alti[i]

    weat_his_amb = DataFrame!(CSV.File("./data/" * site * "/histo_30_amb.csv"))
    weat_his_ele = DataFrame!(CSV.File("./data/" * site * "/histo_30_ele.csv"))
    weat_fut_amb = DataFrame!(CSV.File("./data/" * site * "/rcp85_30_amb.csv"))
    weat_fut_ele = DataFrame!(CSV.File("./data/" * site * "/rcp85_30_ele.csv"))

    opt_his_amb = DataFrame!(CSV.File("./data/" * site * "/opt_his_amb.txt", delim="\t", header=0))
    opt_his_ele = DataFrame!(CSV.File("./data/" * site * "/opt_his_ele.txt", delim="\t", header=0))
    opt_fut_amb = DataFrame!(CSV.File("./data/" * site * "/opt_fut_amb.txt", delim="\t", header=0))
    opt_fut_ele = DataFrame!(CSV.File("./data/" * site * "/opt_fut_ele.txt", delim="\t", header=0))

    year_his = 1971:2005
    year_fut = 2065:2099

    cica_his_amb = GetCiCas(weat_his_amb, year_his, opt_his_amb, lati, long, alti)
    cica_his_ele = GetCiCas(weat_his_ele, year_his, opt_his_ele, lati, long, alti)
    cica_fut_amb = GetCiCas(weat_fut_amb, year_fut, opt_fut_amb, lati, long, alti)
    cica_fut_ele = GetCiCas(weat_fut_ele, year_fut, opt_fut_ele, lati, long, alti)

    writedlm("./data/" * site * "/cicas.txt", [cica_his_amb cica_his_ele cica_fut_amb cica_fut_ele])
end
