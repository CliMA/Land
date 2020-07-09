# include("./scripts/leaf investment/yujie_111_optima_sensitivity_cluster_weather_gs.jl")

list_loca = DataFrame(CSV.File("./data/grid_info.csv"))
list_site = list_loca.grid_id

for i in 1:length(list_site)
    site = list_site[i]
    println(i, " ", site)
    grow_name = homedir() * "/Data/USAPixelSites/GS_" * site * ".csv"
    weat_name = homedir() * "/Data/USAPixelSites/weather_" * site * ".csv"
    grow_data = DataFrame!(CSV.File(grow_name))
    weat_data = DataFrame!(CSV.File(weat_name))

    year_hist = grow_data.Year
    day_start = grow_data.Start_day
    day_end   = grow_data.End_day
    day_gs    = day_end - day_start

    list_loca.MGSL[i] = mean(day_gs)
end

CSV.write("./data/grid_info.csv", list_loca)

#=
list_fixx = CSV.read("./data/0000.csv")
list_loca = CSV.read("./data/grid_info.csv")

for i in 1:length(list_fixx.X)
    id_x = list_fixx.X[i]
    id_y = list_fixx.Y[i]

    # determine where it is
    mask = (list_loca.lat_id .== id_x) .* (list_loca.lon_id .== id_y)
    for nsite in list_loca.indx[mask]
        println(nsite)
        site = list_site[nsite]
        grow_name = homedir() * "/Data/USAPixelSites/GS_" * site * ".csv"
        weat_name = homedir() * "/Data/USAPixelSites/weather_" * site * ".csv"
        grow_data = CSV.read(grow_name)
        weat_data = CSV.read(weat_name)

        year_hist = grow_data.Year
        day_start = grow_data.Start_day
        day_end   = grow_data.End_day
    
        list_gsp = []
        list_gst = []
        for i in 1:length(year_hist)
            year = year_hist[i]
            days = day_start[i]
            daye = day_end[i]
            mask = (weat_data.Year .== year) .* (weat_data.Day .>= days) .* (weat_data.Day .<= daye)
            pric = sum(weat_data.Rain[mask])
            temp = mean(weat_data.Tair[mask])
            append!(list_gsp, pric)
            append!(list_gst, temp)
        end
        mgsp = mean(list_gsp)
        mgst = mean(list_gst)
        println("\tMGSP is ", mgsp)
        println("\tMGST is ", mgst)
        list_loca.MGSP[nsite] = mgsp
        list_loca.MGST[nsite] = mgst

        opti_name = "./data/USA/" * site * ".txt"
        opti_data = CSV.read(opti_name, delim="\t", header=["laba", "vmax", "gcp", "cica"])
        list_loca.LABA[nsite] = mean(opti_data.laba)
        list_loca.VMAX[nsite] = mean(opti_data.vmax)
        list_loca.GCP[nsite]  = mean(opti_data.gcp )
        list_loca.CICA[nsite] = mean(opti_data.cica)
        list_loca.done[nsite] = 1
    end
end
=#

#=
# function to run through sites and years
list_loca = CSV.read("./data/grid_info.csv")
list_site = list_loca.grid_id

#list_loca.MGSP = 0.0
#list_loca.MGST = 0.0
#list_loca.GCP  = 0.0
#list_loca.LABA = 0.0
#list_loca.VMAX = 0.0
#list_loca.CICA = 0.0

for nsite in list_loca.indx[list_loca.done .== 1]
    println(nsite, " in ", length(list_site))
    site = list_site[nsite]
    grow_name = "D:/FIA_Project/WEATHER_DATA/GS_" * site * ".csv"
    weat_name = "D:/FIA_Project/WEATHER_DATA/weather_" * site * ".csv"
    grow_data = CSV.read(grow_name)
    weat_data = CSV.read(weat_name)

    year_hist = grow_data.Year
    day_start = grow_data.Start_day
    day_end   = grow_data.End_day

    list_gsp = []
    list_gst = []
    for i in 1:length(year_hist)
        year = year_hist[i]
        days = day_start[i]
        daye = day_end[i]
        mask = (weat_data.Year .== year) .* (weat_data.Day .>= days) .* (weat_data.Day .<= daye)
        pric = sum(weat_data.Rain[mask])
        temp = mean(weat_data.Tair[mask])
        append!(list_gsp, pric)
        append!(list_gst, temp)
    end
    mgsp = mean(list_gsp)
    mgst = mean(list_gst)
    println("\tPrecipitation is ", mgsp)
    println("\tMean GS Temperature is ", mgst)
    list_loca.MGSP[nsite] = mgsp
    list_loca.MGST[nsite] = mgst
end

CSV.write("./data/grid_info.csv", list_loca)
=#

#=
list_loca = CSV.read("./data/location.csv")
list_file = ["C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/Athens/Athens_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/BoulderMountain/Boulder_Mountain_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/BrasherFalls/Brasher_Falls_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/Durango/Durango_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/Eugene/Eugene_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/Fairbanks/Fairbanks_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/Flagstaff/Flagstaff_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/Hattiesburg/Hattiesburg_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/InternationalFalls/International_Falls_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/Lexington/Lexington_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/McCall/McCall_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/Mio/Mio_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/Missoula/Missoula_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/Moosehorn/Moosehorn_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/Northway/Northway_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/Olustee/Olustee_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/RedFeatherLakes/Red_Feather_Lakes_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/SedroWoolley/Sedro_Woolley_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/Tillamook/Tillamook_histo_sim0_GS.csv",
             "C:/Users/jesin/Documents/MEGA/Research/SPAC - LeafOptimization/6 shifts/Trinity/Trinity_histo_sim0_GS.csv"]

for nsite in 1:20
    println(nsite, " in ", length(list_file))
    weat_name = list_file[nsite]
    weat_data = CSV.read(weat_name)

    year_hist = 1971:2005

    list_gsp = []
    list_gst = []
    for year in year_hist
        mask = weat_data.Year .== year
        pric = sum(weat_data.Rain[mask])
        temp = mean(weat_data.Tair[mask])
        append!(list_gsp, pric)
        append!(list_gst, temp)
    end
    mgsp = mean(list_gsp)
    mgst = mean(list_gst)
    println("\tPrecipitation is ", mgsp)
    println("\tMean GS Temperature is ", mgst)
    list_loca.MGSP[nsite] = mgsp
    list_loca.MGST[nsite] = mgst
end

CSV.write("./data/location.csv", list_loca)
=#