# include("./scripts/leaf investment/yujie_111_optima_sensitivity_cluster_weather_gpp.jl")

list_loca = DataFrame!(CSV.File("./data/grid_info.csv"))
list_site = list_loca.grid_id

@everywhere function get_mean_gpp_thread(param)
    count     = param[1]
    site      = param[2]
    node_opti = Yujie111Init()

    grow_name = homedir() * "/Data/USAPixelSites/GS_"      * site * ".csv"
    weat_name = homedir() * "/Data/USAPixelSites/weather_" * site * ".csv"
    opti_name = "./data/USA/" * site * ".txt"
    grow_data = DataFrame!(CSV.File(grow_name))
    weat_data = DataFrame!(CSV.File(weat_name))
    opti_data = DataFrame!(CSV.File(opti_name, delim="\t", header=["laba", "vmax", "gcp", "cica"]))

    year_hist = grow_data.Year
    day_start = grow_data.Start_day
    day_end   = grow_data.End_day

    list_gpp  = []
    for i in 1:length(year_hist)
        node_temp     = deepcopy(node_opti)
        year          = year_hist[i]
        days          = day_start[i]
        daye          = day_end[i]
        laba          = opti_data.laba[i]
        vmax          = opti_data.vmax[i]
        mask          = (weat_data.Year .== year) .* (weat_data.Day .>= days) .* (weat_data.Day .<= daye)
        weat          = weat_data[mask,:]
        weat[!,:CO2] .= 40.0
        Yujie111UpdateLeaf(node_temp, laba, vmax)
        Yujie111UpdateSoilFromSWC(node_temp, 1.0)
        gpp           = Yujie111GetAnnualGPP(node_temp, weat)
        push!(list_gpp, gpp)
    end

    mean_gpp = mean(list_gpp)
    println("Mean GPP of the " * count * " th data is ", mean_gpp)
    return mean_gpp
end

# generate list for pmap
params = []
for i in 1:length(list_site)
    push!(params, [string(i) * "/" * string(length(list_site)), list_site[i]])
end

results = pmap(get_mean_gpp_thread, params)

list_loca.GPP = results

CSV.write("./data/grid_info.csv", list_loca)