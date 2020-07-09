# include("./scripts/leaf investment/yujie_111_optima_sensitivity_cluster_weather.jl")

# for current ambient
function GenerateOptimaUSA(weather, years, day_s, day_e, filename, d_lati, d_long, d_alti)
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

    @everywhere function GetOptimaForYear(yeardays)
        #println(yeardays)
        year     = yeardays[1]
        dayst    = yeardays[2]
        dayen    = yeardays[3]
        tmp_node = deepcopy(opt_node)
        ini_laba = opt_node.laba
        ini_vmax = opt_node.v_max
        println("Start year ", year, "!")
        weat_mask = (weat_years.Year .== year)
        weat_year = weat_years[weat_mask,:]
        mask_days = (weat_year.Day .>= dayst) .* (weat_year.Day .<= dayen)
        weat_days = weat_year[mask_days,:]
        weat_days[!, :CO2] .= 40.0
        opt_laba,opt_vmax,opt_prof = Yujie111GetOptimalInvestment(tmp_node, weat_days, ini_laba, ini_vmax, max_vmax, false)
        Yujie111UpdateLeaf(opt_node, opt_laba, opt_vmax)
        opt_cica = Yujie111GetAnnualCiCa(opt_node, weat_days)
        println("Finished year ", year, "!")
        return [opt_laba opt_vmax opt_prof opt_cica[1]]
    end

    # run the simulations
    yeardays  = []
    for i in 1:length(years)
        append!(yeardays, [[years[i] day_s[i] day_e[i]]])
    end
    #println(yeardays)
    result = pmap(GetOptimaForYear, yeardays)
    writedlm(filename, result)

    return result
end




# function to run through several sites
list_loca = DataFrame!( CSV.File("./data/grid_info.csv") )
list_site = list_loca.grid_id
list_lati = list_loca.lat
list_long = list_loca.lon
list_alti = list_loca.elev

# generate random sites
total = length(list_loca.done)
udone = sum( list_loca.done .== 0 )
mask  = list_loca.indx[ (list_loca.done .== 0) ]
mode  = "order"
numb  = 500

test_site = []
if udone > numb
    if mode == "random"
        test_site = randsubseq(mask, 1.0*numb/udone)
    else
        temp_loca = list_loca[mask,:]
        test_site = temp_loca[1:numb,:].indx
    end
else
    test_site = mask
end

println("Predicted runs of simulation is ", numb, " out of ", udone)
println("Actaul runs of simulation is ", length(test_site))
println(test_site)

# iterate through the random sites
for nindx in 1:length(test_site)
    nsite = test_site[nindx]
    println("The ", nindx, "th out of ", length(test_site),  " simulations")
    site = list_site[nsite]
    lati = list_lati[nsite]
    long = list_long[nsite]
    alti = list_alti[nsite]

    #grow_name = "C:/Users/jesin/Documents/YUJIE_WEATHER/GS_" * site * ".csv"
    grow_name = homedir() * "/Data/USAPixelSites/GS_" * site * ".csv"
    #weat_name = "C:/Users/jesin/Documents/YUJIE_WEATHER/weather_" * site * ".csv"
    weat_name = homedir() * "/Data/USAPixelSites/weather_" * site * ".csv"
    save_name = "./data/USA/" * site * ".txt"
    println(grow_name)
    println(weat_name)

    grow_data = DataFrame!(CSV.File(grow_name))
    weat_data = DataFrame!(CSV.File(weat_name))

    year_hist = grow_data.Year
    day_start = grow_data.Start_day
    day_end   = grow_data.End_day

    # update the weather data
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

    tmp_simu = GenerateOptimaUSA(weat_name, year_hist, day_start, day_end, save_name, lati, long, alti)
    laba,vmax,gcp,cica = mean(tmp_simu)

    list_loca.done[nsite] = 1
    list_loca.LABA[nsite] = laba
    list_loca.VMAX[nsite] = vmax
    list_loca.GCP[nsite]  = gcp
    list_loca.CICA[nsite] = cica

    # update the list_loca
    CSV.write("./data/grid_info.csv", list_loca)
end
