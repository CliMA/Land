"""
    Tree{FT<:AbstractFloat}

A Tree type which, by default, includes
 - 5 root layers
 - a trunk
 - 20 branches
 - canopy layers

Each canopy layer includes 325 leaves (`36*9` sunlit leaves and 1 shaded leaf).
"""
Base.@kwdef mutable struct Tree{FT<:AbstractFloat}
    # tree information
    "age `[year]`"
    age::Int = 10
    "basal area `[m²]`"
    ba ::FT  = FT( 0.1)
    "ground area `[m²]`"
    ga ::FT  = FT(50.0)
    "tree height `[m]`"
    h  ::FT  = FT( 8.0)
    # tree formation from root to leaves
    "root struct with 5 layers by default"
    roots ::Root   = Root{FT,5}()
    "trunk struct"
    trunk ::Stem   = Stem{FT}()
    "branch struct with 20 layers by default"
    branch::Branch = Branch{FT,20}()
    "canopy struct"
    canopy::Canopy = Canopy{FT,20,325}()
end
