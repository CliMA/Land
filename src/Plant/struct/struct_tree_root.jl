# struct for root layer
Base.@kwdef mutable struct struct_tree_root
    # root structure
    z_lo   ::Float32 = -1.0    # m | root depth
    z_hi   ::Float32 = 0.0     # m | the location for upper root
    r_layer::Float32 = 1.0     # m | average horizental root length
    f_layer::Float32 = 1.0     #   | fraction of biomass and k_max

    # hydraulic parameters
    b     ::Float32 = 2.0       # MPa                  | Weibull B
    c     ::Float32 = 5.0       #                      | Weibull C
    k_max ::Float32 = 1.5625    # mol s^-1 MPa^-1      | maximal hydraulic conductance, = 2700.0 Kg h-1 MPa-1 m^-2 basal area
    k_s   ::Float32 = 15.625    # mol s^-1 MPa^-1 m^-2 | maximal hydraulic conductivity per cross section basal area per root depth
    k_rhiz::Float32 = 1.5e10    # mol s^-1 MPa^-1      | maximal rhizosphere conductance, = 9.72e12 Kg h-1 MPa-1 m^-2 basal area

    # soil parameters
    p_ups    ::Float32        = 0.0                  # MPa | soil water potential (upstream)
    soil_a   ::Float32        = 602.0419             #     | soil texture parameter
    soil_m   ::Float32        = 1.0 - 1.0/1.48       #     | soil texture parameter
    soil_n   ::Float32        = 1.48                 #     | soil texture parameter
    soil_msc ::Float32        = 0.4                  #     | maximal soil volumatric water content, vary with rock fraction
    soil_rwc ::Float32        = 1.0                  #     | soil relative water content
    soil_type::AbstractString = "Sandy Clay Loam"    #     | soil texture
    soil_vwc ::Float32        = 0.4                  #     | soil volumatric water content
    t_soil   ::Float32        = 298.15               # K   | soil temperature

    # flows and pressures (need to be updated with time)
    p_rhiz::Float32 = 0.0    # MPa      | water potential at the root-rhizosphere interface
    p_dos ::Float32 = 0.0    # MPa      | xylem water pressure at the tree base (downstream)
    q     ::Float32 = 0.0    # mol s^-1 | flow rate in the xylem

    # pressure, k, and p_history profile
    k_element::Array{Float32,1} = ones(10) * 83.33     # mol s^-1 MPa^-1 | a list of trunk k_max per element
    p_element::Array{Float32,1} = zeros(10)            # MPa             | a list of trunk xylem pressure per element
    p_history::Array{Float32,1} = zeros(10)            # MPa             | a list of trunk xylem pressure history per element
    t_element::Array{Float32,1} = ones(10) * 298.13    # K               | a list of stem temperature for each element
    z_element::Array{Float32,1} = ones(10) * 0.1       # m               | a list of trunk height per element
end




# function to create root list
function create_root_list(n=5, z_lo=-1.0)
    # create a list of root layers
    root_list = [struct_tree_root() for i in 1:n]

    # iterate through the root list
    for i in 1:n
        # update the k_max for each layer, k_max = k_s * area / delta_h / n
        root_list[i].f_layer   = 1.0 / n
        root_list[i].k_max     = 1.5625 / n
        root_list[i].k_element = ones(10) * 15.625 / n
        # update the root z, z_hi and z_lo for plotting purpose, z-element for computing gravitational pressure drop
        root_list[i].z_hi      = z_lo * (i-1)/n
        root_list[i].z_lo      = z_lo * i/n
        root_list[i].z_element = ones(10) * z_lo * i/n * 0.1
    end

    # return the root list
    return root_list
end




# struct for tree root
Base.@kwdef mutable struct struct_tree_roots
    root_list::Array{struct_tree_root,1} = create_root_list()    # | a list of root layer
end
