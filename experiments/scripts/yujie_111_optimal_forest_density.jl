opt_laba = 1852.0
opt_vmax = 63.7
opt_gaba = 500.0

opt_node = Yujie111Init()
opt_node.gaba = opt_gaba
opt_laba,opt_vmax,opt_prof = Yujie111GetOptimalInvestment(opt_node, weat, opt_laba, opt_vmax)
Yujie111UpdateLeaf(opt_node, opt_laba, opt_vmax)

tmp_node = deepcopy(opt_node)
tmp_prof = Yujie111GetAnnualProfit(tmp_node, weat, 80.0)
max_prof = tmp_prof / tmp_node.gaba
println("Maximal stand profit: ", max_prof)

da = 100.0
while da>1
    global opt_laba,opt_vmax,opt_gaba,opt_node,da,max_prof

    # increase da
    count = 0
    if count<1
        while true
            tmp_node = deepcopy(opt_node)
            tmp_node.gaba = opt_gaba + da
            opt_laba,opt_vmax,opt_prof = Yujie111GetOptimalInvestment(tmp_node, weat, opt_laba, opt_vmax)
            Yujie111UpdateLeaf(tmp_node, opt_laba, opt_vmax)
            Yujie111UpdateSoilFromSWC(node, 1.0)
            tmp_prof  = Yujie111GetAnnualProfit(tmp_node, weat, 80.0)
            tmp_prof /= tmp_node.gaba
            println("Temp stand profit: ", tmp_prof, " @ GaBa: ", tmp_node.gaba)
            if tmp_prof>max_prof
                max_prof  = tmp_prof
                opt_gaba += da
                count    += 1
                println("Maximal stand profit: ", max_prof)
            else
                break
            end
        end
    end

    # decrease da
    if count<1
        while true
            tmp_node = deepcopy(opt_node)
            tmp_node.gaba = opt_gaba - da
            opt_laba,opt_vmax,opt_prof = Yujie111GetOptimalInvestment(tmp_node, weat, opt_laba, opt_vmax)
            Yujie111UpdateLeaf(tmp_node, opt_laba, opt_vmax)
            Yujie111UpdateSoilFromSWC(node, 1.0)
            tmp_prof  = Yujie111GetAnnualProfit(tmp_node, weat, 80.0)
            tmp_prof /= tmp_node.gaba
            println("Temp stand profit: ", tmp_prof, " @ GaBa: ", tmp_node.gaba)
            if tmp_prof>max_prof
                max_prof  = tmp_prof
                opt_gaba -= da
                count    += 1
                println("Maximal stand profit: ", max_prof)
            else
                break
            end
        end
    end

    # make da smaller
    da *= 0.1
end

opt_node.gaba = opt_gaba
Yujie111UpdateLeaf(opt_node, opt_laba, opt_vmax)
