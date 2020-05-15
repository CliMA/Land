# struct for a C3 tree
Base.@kwdef mutable struct StructTree{FT<:AbstractFloat}
    # tree information
    age::Int = 10          # year | age
    ba ::FT  = FT( 0.1)    # m²   | basal area
    ga ::FT  = FT(50.0)    # m²   | ground area
    h  ::FT  = FT( 8.0)    # m    | tree height

    # tree formation from root to leaves
    roots ::StructTreeRoot   = StructTreeRoot{FT,5}()           # | root struct with 5 layers by default
    trunk ::StructTreeStem   = StructTreeStem{FT}()             # | trunk struct
    branch::StructTreeBranch = StructTreeBranch{FT,20}()        # | branch struct with 20 layers by default
    canopy::StructTreeCanopy = StructTreeCanopy{FT,20,325}()    # | canopy struct
end
