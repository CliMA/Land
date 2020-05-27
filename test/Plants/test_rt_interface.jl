using Land.Plant

PLT = Plant

@testset "Plant - RT interface" begin
    FT = Float32

    # create a tree with 10 layers of canopy using Plant module
    tree = PLT.Tree{FT, 5, 10, 325}();

    # create 10-layer canopy structs using CanopyRT module
    canopy_rt, canOpt_rt, canRad_rt, arrayOfLeaves = PLT.initialize_rt_module(n_layer=10, LAI=FT(3.0));

    # test the canopy layers
    @test tree.n_canopy == canopy_rt.nlayers
    @test length(tree.canopy_list[1].leaf_list) == length(canRad_rt.absPAR_sunCab[:,:,1])+1

    # update tree canopy information from CanopyRT
    PLT.update_canopy_from_rt_module!(tree, canopy_rt, canOpt_rt, canRad_rt);
    for canopyi in tree.canopy_list
        @test minimum(canopyi.la_list ) > 0
        @test minimum(canopyi.par_list) > 0
        @test minimum(canopyi.t_list  ) > 200
        @test sum(canopyi.la_list) ≈ canopyi.la
    end
end