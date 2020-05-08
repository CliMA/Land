# struct for stem
Base.@kwdef mutable struct struct_tree_stem
    # trunk structure
    r_layer::Float32 = 1.0    # m | horizental radius of the layer
    z_hi   ::Float32 = 3.0    # m | higher height from the ground
    z_lo   ::Float32 = 0.0    # m | lower height from the ground

    # hydraulic parameter
    b    ::Float32 = 2.0      # MPa                  | Weibull B
    c    ::Float32 = 5.0      #                      | Weibull C
    k_max::Float32 = 8.333    # mol s^-1 MPa^-1      | maximal hydraulic conductance, = 5400.0 * 8/3 Kg h-1 MPa-1 m^-2 basal area
    k_s  ::Float32 = 250.0    # mol s^-1 MPa^-1 m^-2 | maximal hydraulic conductivity per cross section area per tree height

    # flows and pressures
    p_base ::Float32 = 0.0    # MPa      | xylem pressure at the tree basa
    p_joint::Float32 = 0.0    # MPa      | xylem pressure at the trunk-stem or stem-leaf joint
    q      ::Float32 = 0.0    # mol s^-1 | flow rate in the xylem

    # pressure, k, and p_history profile
    k_element::Array{Float32,1} = ones(10) * 83.33     # mol s^-1 MPa^-1 | a list of trunk k_max per element
    p_element::Array{Float32,1} = zeros(10)            # MPa             | a list of trunk xylem pressure per element
    p_history::Array{Float32,1} = zeros(10)            # MPa             | a list of trunk xylem pressure history per element
    t_element::Array{Float32,1} = ones(10) * 298.15    # K               | a list of stem temperature for each element
    z_element::Array{Float32,1} = ones(10) * 0.3       # m               | a list of trunk height per element
end




# struct for tree stem
function create_branch_list(n=20, z_lo=3.0, z_hi=8.0)
    # create a list of branches following the default of stem
    branch_list = [struct_tree_stem() for i in 1:n]

    # iterate through the branch for each canopy layers
    for i in 1:n
        # update the k_max for each layer, k_max = k_s * area / delta_h / n
        branch_list[i].k_max     = 5.0 / n
        branch_list[i].k_element = ones(10) * 50.0 / n
        # update the z for each layer, the z_lo and z_hi are meant for plotting, the z_element is meant for computing pressure drop
        branch_list[i].z_lo      = z_lo + (z_hi-z_lo) * (i-1)/n
        branch_list[i].z_hi      = z_lo + (z_hi-z_lo) * i/n
        branch_list[i].z_element = ones(10) * (z_hi-z_lo) * i/n * 0.1
    end

    # return the branch list
    return branch_list
end




# branch struct
Base.@kwdef mutable struct struct_tree_branch
    # branch structure
    branch_list::Array{struct_tree_stem,1} = create_branch_list()    # | a list of struct_tree_stem
end
