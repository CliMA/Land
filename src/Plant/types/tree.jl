"""
    Tree{FT<:AbstractFloat}

A Tree type which, by default, includes
 - 5 root layers
 - a trunk
 - 20 branches
 - 20 canopy layers

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct Tree{FT<:AbstractFloat}
    # tree information
    "Age `[year]`"
    age::Int = 10
    "Basal area `[m²]`"
    ba ::FT  = FT( 0.1)
    "Ground area `[m²]`"
    ga ::FT  = FT(50.0)
    "Tree height `[m]`"
    h  ::FT  = FT( 8.0)

    # tree formation from root to leaves
    "[`Root`](@ref) which contains 5 [`RootLayer`](@ref) by default"
    roots ::Root   = Root{FT,5}()
    "Trunk using type [`Stem`](@ref)"
    trunk ::Stem   = Stem{FT}()
    "[`Branch`](@ref) which contains 20 [`Stem`](@ref) by default"
    branch::Branch = Branch{FT,20}()
    "[`Canopy`](@ref) which contains 20 [`CanopyLayer`](@ref) by default"
    canopy::Canopy = Canopy{FT,20,325}()
end
