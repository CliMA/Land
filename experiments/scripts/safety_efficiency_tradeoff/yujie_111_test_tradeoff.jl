# include("./scripts/safety_efficiency_tradeoff/yujie_111_test_tradeoff.jl")




# define the cluster function
function GetOptimalBs(weather, years, filename, d_lati, d_long, d_alti)
    # initialize the node
    @everywhere max_vmax      = 100.0
    @everywhere tradeoff_node = Yujie111Init()

    @everywhere tradeoff_node.d_lati = $d_lati
    @everywhere tradeoff_node.d_long = $d_long
    @everywhere tradeoff_node.d_alti = $d_alti
    @everywhere Yujie111UpdateSoilFromSWC(tradeoff_node, 1.0)

    # define the thread function
    @everywhere weat_years = CSV.read($weather)

    @everywhere function GetOptimalBForYear(param)
        year = param[1]
        krat = param[2]
        coun = param[3]
        println(coun)
        println(year)
        println(krat)
        weat_mask = (weat_years.Year .== year)
        weat_year = weat_years[weat_mask,:]

        tradeoff_node = Yujie111Init()
        tradeoff_node.d_lati = $d_lati
        tradeoff_node.d_long = $d_long
        tradeoff_node.d_alti = $d_alti
        Yujie111UpdateSoilFromSWC(tradeoff_node, 1.0)

        tradeoff_node.k_root *= krat
        tradeoff_node.k_stem *= krat
        tradeoff_node.k_sla  *= krat
        Yujie111UpdateLeaf(tradeoff_node, 1000.0, 70.0)
        
        # optimize B, LABA, and VMAX
        last_gcp = -Inf
        tmp_node = deepcopy(tradeoff_node)
        curr_gcp = Yujie111GetAnnualProfit(tmp_node, weat_year)
        while true
            # change LABA
            min_laba = max(1E2,tradeoff_node.laba - 500.0)
            max_laba = min(1E5,tradeoff_node.laba + 500.0)
            mid_laba = tradeoff_node.laba
            while (max_laba-min_laba) > 5.0
                mid_laba = 0.5 * (min_laba+max_laba)
                println("\t\tLABA is currently ", mid_laba)
                tmp_node = deepcopy(tradeoff_node)
                Yujie111UpdateLeaf(tmp_node, mid_laba, tradeoff_node.v_max)
                tmp_gcp0 = Yujie111GetAnnualProfit(tmp_node, weat_year)
                
                tmp_node = deepcopy(tradeoff_node)
                Yujie111UpdateLeaf(tmp_node, mid_laba+5.0, tradeoff_node.v_max)
                tmp_gcp1 = Yujie111GetAnnualProfit(tmp_node, weat_year)
                println("\t\t\tGCP at  lower LABA is ", tmp_gcp0)
                println("\t\t\tGCP at higher LABA is ", tmp_gcp1)

                if tmp_gcp1 > tmp_gcp0
                    min_laba = mid_laba
                else
                    max_laba = mid_laba
                end
            end
            Yujie111UpdateLeaf(tradeoff_node, mid_laba, tradeoff_node.v_max)
            tmp_node = deepcopy(tradeoff_node)
            curr_gcp = Yujie111GetAnnualProfit(tmp_node, weat_year)
            println("\tLast GCP is ", last_gcp)
            println("\tCurr GCP is ", curr_gcp)

            # change VMAX
            min_vmax = max( 0.0001, tradeoff_node.v_max - 20.0)
            max_vmax = min(99.8999, tradeoff_node.v_max + 20.0)
            mid_vmax = tradeoff_node.v_max
            while (max_vmax-min_vmax) > 0.1
                mid_vmax = 0.5 * (min_vmax+max_vmax)
                println("\t\tVcmax is currently ", mid_vmax)
                tmp_node = deepcopy(tradeoff_node)
                Yujie111UpdateLeaf(tmp_node, tradeoff_node.laba, mid_vmax)
                tmp_gcp0 = Yujie111GetAnnualProfit(tmp_node, weat_year)
                
                tmp_node = deepcopy(tradeoff_node)
                Yujie111UpdateLeaf(tmp_node, tradeoff_node.laba, mid_vmax+0.1)
                tmp_gcp1 = Yujie111GetAnnualProfit(tmp_node, weat_year)
                println("\t\t\tGCP at  lower Vcmax is ", tmp_gcp0)
                println("\t\t\tGCP at higher Vcmax is ", tmp_gcp1)

                if tmp_gcp1 > tmp_gcp0
                    min_vmax = mid_vmax
                else
                    max_vmax = mid_vmax
                end
            end
            Yujie111UpdateLeaf(tradeoff_node, tradeoff_node.laba, mid_vmax)
            tmp_node = deepcopy(tradeoff_node)
            curr_gcp = Yujie111GetAnnualProfit(tmp_node, weat_year)
            println("\tLast GCP is ", last_gcp)
            println("\tCurr GCP is ", curr_gcp)

            # change B
            min_b = max( 0.2,tradeoff_node.b_leaf - 0.5)
            max_b = min(10.0,tradeoff_node.b_leaf + 0.5)
            mid_b = tradeoff_node.b_leaf
            while (max_b-min_b) > 0.01
                mid_b = 0.5 * (min_b+max_b)
                println("\t\tB is currently ", mid_b)
                tmp_node = deepcopy(tradeoff_node)
                tmp_node.b_root = mid_b
                tmp_node.b_stem = mid_b
                tmp_node.b_leaf = mid_b
                tmp_gcp0 = Yujie111GetAnnualProfit(tmp_node, weat_year)
                
                tmp_node = deepcopy(tradeoff_node)
                tmp_node.b_root = mid_b + 0.01
                tmp_node.b_stem = mid_b + 0.01
                tmp_node.b_leaf = mid_b + 0.01
                tmp_gcp1 = Yujie111GetAnnualProfit(tmp_node, weat_year)
                println("\t\t\tGCP at  lower B is ", tmp_gcp0)
                println("\t\t\tGCP at higher B is ", tmp_gcp1)
                
                if tmp_gcp1 > tmp_gcp0
                    min_b = mid_b
                else
                    max_b = mid_b
                end
            end
            
            tradeoff_node.b_root = mid_b
            tradeoff_node.b_stem = mid_b
            tradeoff_node.b_leaf = mid_b
            tmp_node = deepcopy(tradeoff_node)
            curr_gcp = Yujie111GetAnnualProfit(tmp_node, weat_year)
            println("\tLast GCP is ", last_gcp)
            println("\tCurr GCP is ", curr_gcp)

            # judge to break
            if (curr_gcp - last_gcp) < 1.0
                break
            else
                last_gcp = curr_gcp
            end
        end
        return [year krat tradeoff_node.b_stem tradeoff_node.laba tradeoff_node.v_max curr_gcp]
    end

    # run simulations
    list_param = []
    count = 0
    for i in 1:length(years)
        for ratio_to_vary in 0.3:0.3:3.0
            count += 1
            push!(list_param, [years[i] ratio_to_vary count])
        end
    end
    result = pmap(GetOptimalBForYear, list_param)
    CSV.write(filename, DataFrame( vcat(result...) ), header=["YEAR","KRAT","OPTB","LABA","VMAX","GCP"])
end




# function to run through several sites
list_loca = CSV.read("./data/location.csv")
list_site = list_loca.Site
list_lati = list_loca.Latitude
list_long = list_loca.Longitude
list_alti = list_loca.Elevation

# Durango      4
# Flagstaff    7 done
# Hattiesburg  8
# Northway    15
# Trinity     20

for i in [4,8,15,20]
    println(i)
    
    site = list_site[i]
    lati = list_lati[i]
    long = list_long[i]
    alti = list_alti[i]

    weat_his_amb =  homedir() * "\\Documents\\MEGA\\Research\\SPAC - LeafOptimization\\shifts\\" * site * "\\" * site *"_histo_sim0_GS.csv"
    weat_years   = CSV.read(weat_his_amb)

    year_his = 1971:2005
    file_to_write = "./data/" * site * "_amb.csv"

    GetOptimalBs(weat_his_amb, year_his, file_to_write, lati, long, alti)
end