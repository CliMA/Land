# Use more advanced algorithm?




function Yujie111GetOptimalInvestment(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            weat::DataFrame,
            ini_laba=2000.0,
            ini_vmax=60.0,
            max_vmax=100.0,
            displaying=true
) where {FT<:AbstractFloat}
    opt_laba = ini_laba
    opt_vmax = ini_vmax
    tmp_node = deepcopy(node)
    Yujie111UpdateLeaf(tmp_node, photo_set, opt_laba, opt_vmax)
    opt_prof = Yujie111GetAnnualProfit(tmp_node, photo_set, weat, max_vmax)
    if displaying
        println("        Optimal LaBa: ", opt_laba)
        println("        Optimal Vcmax: ", opt_vmax)
        println("    Optimal Profit: ", opt_prof)
    end

    # optimize laba and vmax
    df_la = 100.0
    df_vc = 10.0

    list_laba = opt_laba
    list_vmax = opt_vmax
    list_prof = opt_prof

    while (df_la>0.9) && (df_vc>0.09)
        # increase the laba
        count = 0
        if count<1
            while true
                tmp_node = deepcopy(node)
                Yujie111UpdateLeaf(tmp_node, photo_set, opt_laba+df_la, opt_vmax)
                prof = Yujie111GetAnnualProfit(tmp_node, photo_set, weat, max_vmax)
                if prof>opt_prof
                    count    += 1
                    opt_prof  = prof
                    opt_laba += df_la
                    if displaying
                        println("        Optimal LaBa: ", opt_laba)
                        println("    Optimal Profit: ", opt_prof)
                    end
                    list_laba = [list_laba; opt_laba]
                    list_vmax = [list_vmax; opt_vmax]
                    list_prof = [list_prof; opt_prof]
                else
                    break
                end
            end
        end

        # decrease the laba
        if count<1
            while true
                tmp_node = deepcopy(node)
                Yujie111UpdateLeaf(tmp_node, photo_set, max(1.0, opt_laba-df_la), opt_vmax)
                prof = Yujie111GetAnnualProfit(tmp_node, photo_set, weat, max_vmax)
                if prof>opt_prof
                    count    += 1
                    opt_prof  = prof
                    opt_laba -= df_la
                    if opt_laba <= 1.0
                        opt_laba = 1.0
                    end
                    if displaying
                        println("        Optimal LaBa: ", opt_laba)
                        println("    Optimal Profit: ", opt_prof)
                    end
                    list_laba = [list_laba; opt_laba]
                    list_vmax = [list_vmax; opt_vmax]
                    list_prof = [list_prof; opt_prof]
                else
                    break
                end
            end
        end

        # increase the vcmax
        count = 0
        if count<1
            while true
                tmp_node = deepcopy(node)
                Yujie111UpdateLeaf(tmp_node, photo_set, opt_laba, min(opt_vmax+df_vc,max_vmax-1E-6))
                prof = Yujie111GetAnnualProfit(tmp_node, photo_set, weat, max_vmax)
                if prof>opt_prof
                    count    += 1
                    opt_prof  = prof
                    opt_vmax += df_vc
                    if opt_vmax>=max_vmax
                        opt_vmax = max_vmax-1E-6
                    end
                    if displaying
                        println("        Optimal Vcmax: ", opt_vmax)
                        println("    Optimal Profit: ", opt_prof)
                    end
                    list_laba = [list_laba; opt_laba]
                    list_vmax = [list_vmax; opt_vmax]
                    list_prof = [list_prof; opt_prof]
                else
                    break
                end
            end
        end

        # decrease the vmax
        if count<1
            while true
                tmp_node = deepcopy(node)
                Yujie111UpdateLeaf(tmp_node, photo_set, opt_laba, opt_vmax-df_vc)
                prof = Yujie111GetAnnualProfit(tmp_node, photo_set, weat, max_vmax)
                if prof>opt_prof
                    count    += 1
                    opt_prof  = prof
                    opt_vmax -= df_vc
                    if displaying
                        println("        Optimal Vcmax: ", opt_vmax)
                        println("    Optimal Profit: ", opt_prof)
                    end
                    list_laba = [list_laba; opt_laba]
                    list_vmax = [list_vmax; opt_vmax]
                    list_prof = [list_prof; opt_prof]
                else
                    break
                end
            end
        end

        # only if count==0
        if count<1
            df_la *= 0.1
            df_vc *= 0.1
        end
    end

    # return the optima
    return opt_laba,opt_vmax,opt_prof
end
