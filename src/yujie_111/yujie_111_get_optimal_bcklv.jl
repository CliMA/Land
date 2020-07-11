# TODO not yet fixed or tested!
# To be fixed later




function Yujie111GetOptimalBCKLV(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            weat::DataFrame,
            selection,
            max_vmax=100.0,
            displaying=true
) where {FT<:AbstractFloat}
    # determine if optimizing B,C,K,L,V
    yesb,yesc,yesk,yesl,yesv = false,false,false,false,false
    if occursin("B", selection) || occursin("b", selection)
        yesb = true
    end
    if occursin("C", selection) || occursin("c", selection)
        yesc = true
    end
    if occursin("K", selection) || occursin("k", selection)
        yesk = true
    end
    if occursin("L", selection) || occursin("l", selection)
        yesl = true
    end
    if occursin("V", selection) || occursin("v", selection)
        yesv = true
    end

    # optimizing selected traits
    while true
        # use count to break
        count = 0

        # initialize the simulations
        node_temp = deepcopy(node)
        gscp_node = Yujie111GetAnnualProfit(node_temp, weat, max_vmax)

        # optimize V
        if yesv
            diff_v = 10.0
            if displaying
                println("\tdiff_v   returns to ", diff_v)
            end
            while diff_v>0.09
                # increase vmax
                while true
                    node_temp = deepcopy(node)
                    Yujie111UpdateLeaf(node_temp, photo_set, node_temp.laba, min(node_temp.ps.Vcmax25+diff_v, max_vmax-1E-6))
                    gscp_incr = Yujie111GetAnnualProfit(node_temp, weat, max_vmax)
                    if gscp_incr>gscp_node
                        gscp_node = gscp_incr
                        Yujie111UpdateLeaf(node, photo_set, node.laba, min(node.ps.Vcmax25+diff_v, max_vmax-1E-6))
                        if displaying
                            println("\t\tOptimal Vcmax increases by ", diff_v, " to ", node.ps.Vcmax25)
                        end
                        count += 1
                    else
                        break
                    end
                end
                # decrease vmax
                while true
                    node_temp = deepcopy(node)
                    Yujie111UpdateLeaf(node_temp, photo_set, node_temp.laba, max(node_temp.ps.Vcmax25-diff_v, 0.1))
                    gscp_decr = Yujie111GetAnnualProfit(node_temp, weat, max_vmax)
                    if gscp_decr>gscp_node
                        gscp_node = gscp_decr
                        Yujie111UpdateLeaf(node, photo_set, node.laba, max(node.ps.Vcmax25-diff_v, 0.1))
                        if displaying
                            println("\t\tOptimal Vcmax decreases by ", diff_v, " to ", node.ps.Vcmax25)
                        end
                        count += 1
                    else
                        break
                    end
                end
                # decrease diff_v
                diff_v *= 0.1
                if displaying
                    println("\tdiff_v decreases to ", diff_v)
                end
            end
        end

        # optimize L
        if yesl
            diff_l = 100.0
            if displaying
                println("\tdiff_l   returns to ", diff_l)
            end
            while diff_l>0.9
                # increase laba
                while true
                    node_temp = deepcopy(node)
                    Yujie111UpdateLeaf(node_temp, photo_set, node_temp.laba+diff_l, node_temp.ps.Vcmax25)
                    gscp_incr = Yujie111GetAnnualProfit(node_temp, weat, max_vmax)
                    if gscp_incr>gscp_node
                        gscp_node = gscp_incr
                        Yujie111UpdateLeaf(node, photo_set, node.laba+diff_l, node.ps.Vcmax25)
                        if displaying
                            println("\t\tOptimal LA:BA increases by ", diff_l, " to ", node.laba)
                        end
                        count += 1
                    else
                        break
                    end
                end
                # decrease laba
                while true
                    node_temp = deepcopy(node)
                    Yujie111UpdateLeaf(node_temp, photo_set, max(0.1,node_temp.laba-diff_l), node_temp.ps.Vcmax25)
                    gscp_decr = Yujie111GetAnnualProfit(node_temp, weat, max_vmax)
                    if gscp_decr>gscp_node
                        gscp_node = gscp_decr
                        Yujie111UpdateLeaf(node, photo_set, max(0.1,node.laba-diff_l), node.ps.Vcmax25)
                        if displaying
                            println("\t\tOptimal LA:BA decreases by ", diff_l, " to ", node.laba)
                        end
                        count += 1
                    else
                        break
                    end
                end
                # decrease diff_l when
                diff_l *= 0.1
                if displaying
                    println("\tdiff_l decreases to ", diff_l)
                end
            end
        end

        # optimize K
        if yesk
            diff_k = 0.1
            if displaying
                println("\tdiff_k   returns to ", diff_k)
            end
            while diff_k>9E-4
                # increase k
                while true
                    node_temp = deepcopy(node)
                    if node_temp.k_root * (1.0 + diff_k) > 1E5
                        break
                    end
                    node_temp.k_root *= 1.0 + diff_k
                    node_temp.k_stem *= 1.0 + diff_k
                    node_temp.k_sla  *= 1.0 + diff_k
                    Yujie111UpdateLeaf(node_temp, photo_set, node_temp.laba, node_temp.ps.Vcmax25)
                    gscp_incr = Yujie111GetAnnualProfit(node_temp, weat, max_vmax)
                    if gscp_incr>gscp_node
                        gscp_node = gscp_incr
                        node.k_root *= 1.0 + diff_k
                        node.k_stem *= 1.0 + diff_k
                        node.k_sla  *= 1.0 + diff_k
                        Yujie111UpdateLeaf(node, photo_set, node.laba, node.ps.Vcmax25)
                        if displaying
                            println("\t\tOptimal K increases by ", diff_k*100.0, "% to ", node.k_root, " (root)")
                        end
                        count += 1
                    else
                        break
                    end
                end
                # decrease k
                while true
                    node_temp = deepcopy(node)
                    if node_temp.k_root * (1.0 - diff_k) < 10.0
                        break
                    end
                    node_temp.k_root *= 1.0 - diff_k
                    node_temp.k_stem *= 1.0 - diff_k
                    node_temp.k_sla  *= 1.0 - diff_k
                    Yujie111UpdateLeaf(node_temp, photo_set, node_temp.laba, node_temp.ps.Vcmax25)
                    gscp_decr = Yujie111GetAnnualProfit(node_temp, weat, max_vmax)
                    if gscp_decr>gscp_node
                        gscp_node = gscp_decr
                        node.k_root *= 1.0 - diff_k
                        node.k_stem *= 1.0 - diff_k
                        node.k_sla  *= 1.0 - diff_k
                        Yujie111UpdateLeaf(node, photo_set, node.laba, node.ps.Vcmax25)
                        if displaying
                            println("\t\tOptimal K decreases by ", diff_k*100.0, "% to ", node.k_root, " (root)")
                        end
                        count += 1
                    else
                        break
                    end
                end
                # change diff_k
                diff_k *= 0.1
                if displaying
                    println("\tdiff_k decreases to ", diff_k)
                end
            end
        end

        # optimize B
        if yesb
            diff_b = 1.0
            if displaying
                println("\tdiff_b   returns to ", diff_b)
            end
            while diff_b>0.009
                # increase b
                while true
                    node_temp = deepcopy(node)
                    node_temp.b_root += diff_b
                    node_temp.b_stem += diff_b
                    node_temp.b_leaf += diff_b
                    gscp_incr = Yujie111GetAnnualProfit(node_temp, weat, max_vmax)
                    if gscp_incr>gscp_node
                        gscp_node = gscp_incr
                        node.vc_root.b += diff_b
                        node.vc_stem.b += diff_b
                        node.vc_leaf.b += diff_b
                        if displaying
                            println("\t\tOptimal B increases by ", diff_b, " to ", node.vc_root.b)
                        end
                        count += 1
                    else
                        break
                    end
                end
                # decrease b
                while true
                    node_temp = deepcopy(node)
                    node_temp.b_root = max(0.1,node_temp.b_root-diff_b)
                    node_temp.b_stem = max(0.1,node_temp.b_stem-diff_b)
                    node_temp.b_leaf = max(0.1,node_temp.b_leaf-diff_b)
                    gscp_decr = Yujie111GetAnnualProfit(node_temp, weat, max_vmax)
                    if gscp_decr>gscp_node
                        gscp_node = gscp_decr
                        node.vc_root.b = max(0.1,node.vc_root.b-diff_b)
                        node.vc_stem.b = max(0.1,node.vc_stem.b-diff_b)
                        node.vc_leaf.b = max(0.1,node.vc_leaf.b-diff_b)
                        if displaying
                            println("\t\tOptimal B decreases by ", diff_b, " to ", node.vc_root.b)
                        end
                        count += 1
                    else
                        break
                    end
                end
                # decrease diff_b
                diff_b *= 0.1
                if displaying
                    println("\tdiff_b decreases to ", diff_b)
                end
            end
        end

        # this part has not yet been tested because it runs into some serious problem of C<1!
        # optimize C
        if yesc
            diff_c = 1.0
            if displaying
                println("\tdiff_c   returns to ", diff_c)
            end
            while diff_c>0.009
                # increase c
                while true
                    node_temp = deepcopy(node)
                    node_temp.c_root += diff_c
                    node_temp.c_stem += diff_c
                    node_temp.c_leaf += diff_c
                    gscp_incr = Yujie111GetAnnualProfit(node_temp, weat, max_vmax)
                    if gscp_incr>gscp_node
                        gscp_node = gscp_incr
                        node.vc_root.c += diff_c
                        node.vc_stem.c += diff_c
                        node.vc_leaf.c += diff_c
                        if displaying
                            println("\t\tOptimal C increases by ", diff_c, " to ", node.vc_root.c)
                        end
                        count += 1
                    else
                        break
                    end
                end
                # decrease c
                while true
                    node_temp = deepcopy(node)
                    node_temp.c_root = max(1.1,node_temp.c_root-diff_c)
                    node_temp.c_stem = max(1.1,node_temp.c_stem-diff_c)
                    node_temp.c_leaf = max(1.1,node_temp.c_leaf-diff_c)
                    gscp_decr = Yujie111GetAnnualProfit(node_temp, weat, max_vmax)
                    if gscp_decr>gscp_node
                        gscp_node = gscp_decr
                        node.vc_root.c = max(1.1,node.vc_root.c - diff_c)
                        node.vc_stem.c = max(1.1,node.vc_stem.c - diff_c)
                        node.vc_leaf.c = max(1.1,node.vc_leaf.c - diff_c)
                        if displaying
                            println("\t\tOptimal C decreases by ", diff_c, " to ", node.vc_root.c)
                        end
                        count += 1
                    else
                        break
                    end
                end
                # decrease diff_c
                diff_c *= 0.1
                if displaying
                    println("\tdiff_c decreases to ", diff_c)
                end
            end
        end

        # judge whether to break
        if count==0
            break
        end
    end

    # return the results
    return [node.vc_root.b node.vc_root.c node.k_root node.k_stem node.k_sla node.laba node.ps.Vcmax25]
end