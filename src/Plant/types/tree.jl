"""
    Tree{FT<:AbstractFloat}

A Tree type which includes
- `root_no` number of root layers
- a trunk
- `canopy_no` branches
- `leaf_no` canopy layers

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct Tree{FT<:AbstractFloat, root_no, canopy_no, leaf_no}
    # plant information
    "Age `[year]`"
    age::Int = 10
    "Basal area `[m²]`"
    ba ::FT  = FT( 0.1)
    "Ground area `[m²]`"
    ga ::FT  = FT(50.0)
    "Tree height `[m]`"
    h  ::FT  = FT( 8.0)
    "Photosynthesis type C3/C4/CAM"
    photo_type::String = "C3"
    "Photosynthesis model parameter set"
    photo_para_set::C3VcVpJBernacchi = C3VcVpJBernacchi{FT}()
    "Optimization model option"
    stomata_scheme::AbstractStomatalModel = OSMWang()

    # structure information

    # tree formation of root
    "Number of root layers"
    n_root   ::Int                    = root_no
    "List of root fraction"
    root_frac::Array{FT,1}            = ones(FT,5) .* FT(0.2)
    "List of root z lower boundary"
    root_zs  ::Array{FT,1}            = FT.([-0.2,-0.4,-0.6,-0.8,-1.0])
    "List of [`RootLayer`](@ref)"
    root_list::Array{RootLayer{FT},1} = create_root_list(FT,n_root)

    # tree information of trunk and branch
    "Number of canopy layers"
    n_canopy   ::Int               = canopy_no
    "Trunk using type [`Stem`](@ref)"
    trunk      ::Stem              = Stem{FT}()
    "List of [`Stem`](@ref)"
    branch_list::Array{Stem{FT},1} = create_branch_list(FT,n_canopy)

    # tree information of canopy
    "Number of leaves in each canopy layer"
    n_leaf     ::Int                      = leaf_no
    "List of [`CanopyLayer`](@ref)"
    canopy_list::Array{CanopyLayer{FT},1} = [CanopyLayer{FT,n_leaf}() for i in 1:n_canopy]
end
