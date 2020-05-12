# struct for a C3 tree
Base.@kwdef mutable struct StructTree
    # tree information
    age::Int32   = 10      # year | age
    ba ::Float32 = 0.1     # m²   | basal area
    ga ::Float32 = 50.0    # m²   | ground area
    h  ::Float32 = 8.0     # m    | tree height

    # tree formation from root to leaves
    roots ::StructTreeRoot   = StructTreeRoot()      # | root struct with 5 layers by default
    trunk ::StructTreeStem   = StructTreeStem()      # | trunk struct
    branch::StructTreeBranch = StructTreeBranch()    # | branch struct with 20 layers by default
    canopy::StructTreeCanopy = StructTreeCanopy()    # | canopy struct
end
