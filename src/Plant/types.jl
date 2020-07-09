# the Node111 struct means one root, stem, and leaf layer
Base.@kwdef mutable struct PlantSimple{FT<:AbstractFloat}
    # Plant hydraulics system
    hs::TreeSimple{FT} = TreeSimple{FT}();

    # Leaf photosynthesis system
    ps::Leaf{FT} = Leaf{FT}()

    # Surrounding AirLayer
    envir::AirLayer{FT} = AirLayer{FT}()

    # leaf related
    laba ::FT = 1500    # leaf area per basal area
    g_sla::FT = 0.8     # g per leaf area at 25 C
    gaba ::FT = 500     # ground area per basal area
    width::FT = 0.05    # leaf width in m

    # soil water related
    c_curr::FT = 0.476    # theta at current scenario
    h_soil::FT = 2     # soil depth, which is 2X root depth here

    # pressure distributions
    p_soil::FT = 0.0    # soil P
    p_l_su::FT = 0.0    # P at end of sunlit leaf
    p_l_sh::FT = 0.0    # P at end of shade leaf

    # leaf invest related
    c_cons::FT = 2       # leaf construction cost
    c_vmax::FT = 0.02    # leaf cost in photosynthesis capacity

    # environment related
    d_lati::FT = 30     # latitude in degree
    d_long::FT = 116    # longitude in degree
    d_alti::FT = 0      # altitude in m
end
