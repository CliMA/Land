#=
# a customized function to convert number to formatted number
function number_to_string(number::Number; digit::Int=6)
    if digit==0
        return @sprintf("%0.0f", number)
    elseif digit==1
        return @sprintf("%0.1f", number)
    elseif digit==2
        return @sprintf("%0.2f", number)
    elseif digit==3
        return @sprintf("%0.3f", number)
    elseif digit==4
        return @sprintf("%0.4f", number)
    elseif digit==5
        return @sprintf("%0.5f", number)
    else
        return @sprintf("%0.6f", number)
    end
end




# there is an option to save the figure
function visualize_struct_tree(tree::Tree, save::Bool=false)
    # clear plot first
    figure("Tree Visualization",figsize=(16,9), dpi=120)
    clf()

    # display table header
    text(-1.00, tree.h+0.05, L"$P_\mathrm{water}$"            , color=:red   , ha="right", va="bottom", fontsize=8)
    text( 1.00, tree.h+0.05, L"$Q_\mathrm{water}$"            , color=:blue  , ha="left" , va="bottom", fontsize=8)
    text( 1.40, tree.h+0.05, L"$\overline{\mathrm{PAR_{sl}}}$", color=:orange, ha="left" , va="bottom", fontsize=8)
    text( 1.60, tree.h+0.05, L"$\mathrm{PAR_{sh}}$"           , color=:orange, ha="left" , va="bottom", fontsize=8)

    # plot the root system
    for i in 1:length(tree.roots.root_list)
        rooti = tree.roots.root_list[i]

        lx = [         0,          0]
        ly = [rooti.z_hi, rooti.z_lo]
        plot(lx, ly, color=:brown)

        lx = [-rooti.r_layer,          0, rooti.r_layer]
        ly = [ rooti.z_lo   , rooti.z_hi, rooti.z_lo   ]
        plot(lx, ly, color=:brown)

        text(-rooti.r_layer, rooti.z_lo, number_to_string(rooti.p_ups , digit=3), color=:red , ha="right", va="top"   , fontsize=8)
        text(-rooti.r_layer, rooti.z_lo, number_to_string(rooti.p_rhiz, digit=3), color=:red , ha="left" , va="bottom", fontsize=8)
        text( rooti.r_layer, rooti.z_lo, number_to_string(rooti.q     , digit=6), color=:blue, ha="left" , va="top"   , fontsize=8)
    end

    # plot the trunk
    lx = [0, 0, 0]
    ly = [0, tree.trunk.z_lo, tree.trunk.z_hi]
    plot(lx, ly, color=:black)
    text(0, tree.trunk.z_lo    , number_to_string(tree.trunk.p_ups, digit=3), color=:red , ha="right", va="bottom", fontsize=8)
    text(0, tree.trunk.z_hi    , number_to_string(tree.trunk.p_dos, digit=3), color=:red , ha="right", va="top"   , fontsize=8)
    text(0, tree.trunk.z_hi*0.5, number_to_string(tree.trunk.q    , digit=6), color=:blue, ha="left" , va="top"   , fontsize=8)

    # plot the canopy stem and leaf system
    for i in 1:length(tree.branch.branch_list)
        branchi = tree.branch.branch_list[i]
        canopyi = tree.canopy.canopy_list[i]

        # plot the branch vertically
        lx = [           0,            0]
        ly = [branchi.z_lo, branchi.z_hi]
        plot(lx, ly, color=:green)

        # plot the branch horizentally
        lx = [-branchi.r_layer,            0, branchi.r_layer]
        ly = [ branchi.z_hi   , branchi.z_lo, branchi.z_hi   ]
        plot(lx, ly, color=:green)

        # label the stem-leaf joint pressure and flow rate through each branch
        text(-branchi.r_layer, branchi.z_hi, number_to_string(branchi.p_dos, digit=3), color=:red , ha="right", va="top", fontsize=8)
        text( branchi.r_layer, branchi.z_hi, number_to_string(branchi.q    , digit=6), color=:blue, ha="left" , va="top", fontsize=8)

        # calculate and display mean PAR for sunlit leaves and shaded leaves
        par_sl = sum(canopyi.par_list[1:end-1] .* canopyi.la_list[1:end-1]) / sum(canopyi.la_list[1:end-1])
        par_sh = canopyi.par_list[end]
        text(1.40, branchi.z_hi, number_to_string(par_sl, digit=1), color=:orange, ha="left", va="top", fontsize=8)
        text(1.60, branchi.z_hi, number_to_string(par_sh, digit=1), color=:orange, ha="left", va="top", fontsize=8)
    end

    # plot the ground
    plot([-1.2,4], [0,0], "k-", color=:black, linewidth=0.2)

    # set the title and label
    title("Tree system")
    xlabel("Distance from the tree base (m)", fontsize=16)
    ylabel("z (m)"                          , fontsize=16)

    # save the figure
    if save
        savefig("zzz_tree.pdf", bbox_inches="tight")
    end
end
=#