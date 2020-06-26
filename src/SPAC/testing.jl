



function test_diurnal_cycle(Δt::FT) where {FT<:AbstractFloat}
    # define soil and air layer boundaries
    soil_bounds = FT[0, -0.5, -1, -1.5, -2, -3, -4, -5, -6, -7, -8, -9, -10];
    air_bounds  = collect(FT, 0:1:120);

    # create a tree with root-trunk-canopy z of -2, 5, and 10 m
    tree::Tree{FT} = create_tree(FT(-2), FT(5), FT(10), soil_bounds, air_bounds);

    # create canopy accordingly
    angles,arrayOfLeaves,canopy_rt,canOpt_rt,canRad_rt,soil,sunRad_rt,wl_set = 
                initialize_rt_module(n_layer=tree.n_canopy, LAI=FT(3));

    #
end
