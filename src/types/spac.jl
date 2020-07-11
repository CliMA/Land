# the Node111 struct means one root, stem, and leaf layer
Base.@kwdef mutable struct SPACSimple{FT<:AbstractFloat}
    # hydraulic systems for SPACSimple
    hs::TreeSimple{FT} = TreeSimple{FT}()
    "Critical flow rate"
    ec::FT = FT(50);

    # photosynthesis parameters
    ps::Leaf = Leaf{FT}()

    # Surrounding AirLayer
    envir::AirLayer{FT} = AirLayer{FT}()

    # local container for returned results
    container1L::SPACContainer1L{FT} = SPACContainer1L{FT}()
    container2L::SPACContainer2L{FT} = SPACContainer2L{FT}()

    # cica information
    cica_all::FT = FT(0)
    cica_fav::FT = FT(0)

    # opt fs information
    opt_f_sl::FT = FT(0)
    opt_f_sh::FT = FT(0)

    # leaf related
    lai  ::FT = FT(3)   # leaf area index
    laba ::FT = 1500    # leaf area per basal area
    g_sla::FT = 0.8     # g per leaf area at 25 C
    gaba ::FT = 500     # ground area per basal area
    width::FT = 0.05    # leaf width in m

    # soil related
    c_ssat::FT = 0.476    # theta at saturation
    c_curr::FT = 0.476    # theta at current scenario
    k_ssat::FT = 0.009    # kmax at saturation
    p_soil::FT = 0.0    # soil P
    h_soil::FT = 2     # soil depth, which is 2X root depth here

    # leaf invest related
    c_cons::FT = 2       # leaf construction cost
    c_vmax::FT = 0.02    # leaf cost in photosynthesis capacity

    # geography related
    d_lati::FT = 30     # latitude in degree
    d_long::FT = 116    # longitude in degree
    d_alti::FT = 0      # altitude in m
end
