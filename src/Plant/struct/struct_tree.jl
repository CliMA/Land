# struct for a C3 tree
Base.@kwdef mutable struct Tree{FT<:AbstractFloat}
    # tree information
    age::Int = 10          # year | age
    ba ::FT  = FT( 0.1)    # m²   | basal area
    ga ::FT  = FT(50.0)    # m²   | ground area
    h  ::FT  = FT( 8.0)    # m    | tree height

    # tree formation from root to leaves
    roots ::Root   = Root{FT,5}()           # | root struct with 5 layers by default
    trunk ::Stem   = Stem{FT}()             # | trunk struct
    branch::Branch = Branch{FT,20}()        # | branch struct with 20 layers by default
    canopy::Canopy = Canopy{FT,20,325}()    # | canopy struct
end
