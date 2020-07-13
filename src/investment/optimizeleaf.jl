# Use more advanced algorithm?


function optimize_leaf!(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            weat::Array{FT,2};
            displaying=true
) where {FT<:AbstractFloat}
    # 1. use the opt_laba and opt_vmax to initialize
    opt_laba = node.opt_laba;
    opt_vmax = node.opt_vmax;
    tmp_node = deepcopy(node);
    leaf_allocation!(tmp_node, photo_set, opt_laba, opt_vmax);
    opt_prof = annual_profit(tmp_node, photo_set, weat);
    if displaying
        println("    OPT_LABA -> ", opt_laba);
        println("    OPT_VMAX -> ", opt_vmax);
        println("    MAX_GSCP -> ", opt_prof);
    end

    # 2. define the steps for LABA and Vcmax
    step_laba  = FT(100);
    step_vmax  = FT(10);

    # 3. stepwisely optimize the parameters
    while (step_laba > 0.9) && (step_vmax>0.09)
        # 3.1 increase the laba
        count_laba = 0;
        while true
            tmp_node = deepcopy(node);
            tmp_laba = opt_laba + step_laba;
            leaf_allocation!(tmp_node, photo_set, tmp_laba, opt_vmax);
            tmp_prof = annual_profit(tmp_node, photo_set, weat);
            if tmp_prof > opt_prof
                count_laba += 1;
                opt_prof    = tmp_prof;
                opt_laba    = tmp_laba;
                if displaying
                    println("    OPT_LABA -> ", opt_laba)
                    println("    MAX_GSCP -> ", opt_prof)
                end
            else
                break
            end
        end

        # 3.2 decrease the laba
        if count_laba==0
            while true
                tmp_node = deepcopy(node);
                tmp_laba = opt_laba - step_laba;
                tmp_laba <= 0 ? break : nothing;
                leaf_allocation!(tmp_node, photo_set, tmp_laba, opt_vmax);
                tmp_prof = annual_profit(tmp_node, photo_set, weat);
                if tmp_prof > opt_prof
                    count_laba += 1;
                    opt_prof    = tmp_prof;
                    opt_laba    = tmp_laba;
                    if displaying
                        println("    OPT_LABA -> ", opt_laba)
                        println("    MAX_GSCP -> ", opt_prof)
                    end
                else
                    break
                end
            end
        end

        # 3.3 increase the Vcmax25
        count_vmax = 0;
        while true
            tmp_node = deepcopy(node);
            tmp_vmax = opt_vmax + step_vmax;
            tmp_vmax >= tmp_node.maxv ? break : nothing;
            leaf_allocation!(tmp_node, photo_set, opt_laba, tmp_vmax);
            tmp_prof = annual_profit(tmp_node, photo_set, weat);
            if tmp_prof > opt_prof
                count_vmax += 1;
                opt_prof    = tmp_prof;
                opt_vmax    = tmp_vmax;
                if displaying
                    println("    OPT_VMAX -> ", opt_vmax)
                    println("    MAX_GSCP -> ", opt_prof)
                end
            else
                break
            end
        end

        # 3.4 decrease the Vcmax25
        if count_vmax==0
            while true
                tmp_node = deepcopy(node);
                tmp_vmax = opt_vmax - step_vmax;
                tmp_vmax <= 0 ? break : nothing;
                leaf_allocation!(tmp_node, photo_set, opt_laba, tmp_vmax);
                tmp_prof = annual_profit(tmp_node, photo_set, weat);
                if tmp_prof > opt_prof
                    count_vmax += 1;
                    opt_prof    = tmp_prof;
                    opt_vmax    = tmp_vmax;
                    if displaying
                        println("    OPT_VMAX -> ", opt_vmax)
                        println("    MAX_GSCP -> ", opt_prof)
                    end
                else
                    break
                end
            end
        end

        # 3.5 reduce the steps
        if (count_laba==0) && (count_vmax==0)
            step_laba /= 10;
            step_vmax /= 10;
        end
    end

    # 4. update the opt_laba and opt_vmax to node
    leaf_allocation!(node, photo_set, opt_laba, opt_vmax);
    node.opt_laba = opt_laba;
    node.opt_vmax = opt_vmax;

    return nothing
end
