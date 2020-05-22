"""
    Tree{FT<:AbstractFloat}

A Tree type which includes
- `n_root` number of root layers
- a trunk
- `n_canopy` branches
- `n_canopy` canopy layers

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct Tree{FT<:AbstractFloat, n_root, n_canopy, n_leaf}
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

    # tree formation from root to leaves
    "[`Root`](@ref) which contains n_root [`RootLayer`](@ref)"
    roots ::Root   = Root{FT,n_root}()
    "Trunk using type [`Stem`](@ref)"
    trunk ::Stem   = Stem{FT}()
    "[`Branch`](@ref) which contains n_canopy [`Stem`](@ref)"
    branch::Branch = Branch{FT,n_canopy}()
    "[`Canopy`](@ref) which contains n_canopy [`CanopyLayer`](@ref)"
    canopy::Canopy = Canopy{FT,n_canopy,n_leaf}()
end
