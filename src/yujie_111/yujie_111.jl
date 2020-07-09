# the Node111 struct means one root, stem, and leaf layer
Base.@kwdef mutable struct Yujie111{FT<:AbstractFloat}
    # hydraulic systems for Yujie111
    root = RootHydraulics{FT}();
    stem = StemHydraulics{FT}();
    leaf = LeafHydraulics{FT}();

    # hydraulic conductance for tree and each element
    k_tree::FT = 2500     # tree maximal K
    k_rhiz::FT = 5e12     # rhizosphere maximal K
    k_root::FT = 5000     # root maximal K
    k_stem::FT = 10000    # stem maximal K
    k_leaf::FT = 10000    # leaf maximal K

    # leaf area related
    k_sla::FT = 2       # leaf maximal K per leaf area
    laba ::FT = 5000    # leaf area per basal area
    g_sla::FT = 0.8     # g per leaf area at 25 C
    gmax ::FT = laba * g_sla * 3600.0 * 18.0 * 0.001   # maximal gmax per tree at 25 C
    gaba ::FT = 500     # ground area per basal area
    width::FT = 0.05    # leaf width in m

    # vulnerability parameters
    p_ssat::FT = 0.63*9.8*998.0*1E-6    # phi at saturation
    c_ssat::FT = 0.476    # theta at saturation
    c_curr::FT = 0.476    # theta at current scenario
    b_ssat::FT = 8.52     # b value for ths soil
    k_ssat::FT = 0.009    # kmax at saturation

    # vulnerability curves for xylem
    vc_root::WeibullSingle{FT} = WeibullSingle{FT}()
    vc_stem::WeibullSingle{FT} = WeibullSingle{FT}()
    vc_leaf::WeibullSingle{FT} = WeibullSingle{FT}()

    # pressure distributions
    p_soil::FT = 0.0    # soil P
    p_rhiz::FT = 0.0    # P at end of rhizosphere
    p_root::FT = 0.0    # P at end of root
    p_stem::FT = 0.0    # P at end of stem
    p_l_su::FT = 0.0    # P at end of sunlit leaf
    p_l_sh::FT = 0.0    # P at end of shade leaf

    # legacies for root, stem, and leaf
    l_root::Array{FT,2}  = [zeros(20) ones(20)]    # K legacy for root
    l_stem::Array{FT,2}  = [zeros(20) ones(20)]    # K legacy for stem
    l_leaf::Array{FT,2}  = [zeros(20) ones(20)]    # K legacy for leaf

    # photosynthesis parameters
    ps::Leaf = Leaf{FT}()

    # height related
    h_soil::FT = 2     # soil depth, which is 2X root depth here
    h_root::FT = 1     # mean root depth
    h_stem::FT = 10    # stem height
    h_leaf::FT = 0    # leaf heihgt, neglectable in most trees

    # leaf invest related
    c_cons::FT = 2       # leaf construction cost
    c_vmax::FT = 0.02    # leaf cost in photosynthesis capacity

    # environment related
    d_lati::FT = 30     # latitude in degree
    d_long::FT = 116    # longitude in degree
    d_alti::FT = 0      # altitude in m
end
