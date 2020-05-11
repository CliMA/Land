# struct for a C3 tree
Base.@kwdef mutable struct struct_tree
    # tree information
    age::Int32   = 10      # year | age
    ba ::Float32 = 0.1     # m^2  | basal area
    ga ::Float32 = 50.0    # m^2  | ground area
    h  ::Float32 = 8.0     # m    | tree height

    # tree formation from root to leaves
    roots ::struct_tree_roots   = struct_tree_roots()     # | root struct with 5 layers by default
    trunk ::struct_tree_stem    = struct_tree_stem()      # | trunk struct
    branch::struct_tree_branch  = struct_tree_branch()    # | branch struct with 20 layers by default
    canopy ::struct_tree_canopy = struct_tree_canopy()    # | canopy struct
end
