# struct for root layer
Base.@kwdef mutable struct RootLayer{FT<:AbstractFloat}
    # root structure
    z_lo   ::FT = FT(-1.0)    # m | root depth
    z_hi   ::FT = FT( 0.0)    # m | the location for upper root
    r_layer::FT = FT( 1.0)    # m | average horizental root length
    f_layer::FT = FT( 1.0)    #   | fraction of biomass and k_max

    # hydraulic parameters
    b     ::FT = FT( 2.0   )    # MPa               | Weibull B
    c     ::FT = FT( 5.0   )    #                   | Weibull C
    k_max ::FT = FT(1.5625 )    # mol s⁻¹ MPa⁻¹     | maximal hydraulic conductance, = 2700.0 Kg h⁻¹ MPa⁻¹ m⁻² basal area
    k_s   ::FT = FT(15.625 )    # mol s⁻¹ MPa⁻¹ m⁻² | maximal hydraulic conductivity per cross section basal area per root depth
    k_rhiz::FT = FT( 1.5e10)    # mol s⁻¹ MPa⁻¹     | maximal rhizosphere conductance, = 9.72e12 Kg h⁻¹ MPa⁻¹ mm⁻² basal area

    # soil parameters
    p_ups    ::FT     = FT(  0.0   )         # MPa | soil water potential (upstream)
    soil_a   ::FT     = FT(602.0419)         #     | soil texture parameter
    soil_m   ::FT     = FT(  3.2432)         #     | soil texture parameter, 1- 1/n
    soil_n   ::FT     = FT(  1.48  )         #     | soil texture parameter
    soil_msc ::FT     = FT(  0.4   )         #     | maximal soil volumatric water content, vary with rock fraction
    soil_rwc ::FT     = FT(  1.0   )         #     | soil relative water content
    soil_vwc ::FT     = FT(  0.4   )         #     | soil volumatric water content
    t_soil   ::FT     = K_25                 # K   | soil temperature
    soil_type::String = "Sandy Clay Loam"    #     | soil texture

    # flows and pressures (need to be updated with time)
    p_rhiz::FT = FT(0.0)    # MPa     | water potential at the root-rhizosphere interface
    p_dos ::FT = FT(0.0)    # MPa     | xylem water pressure at the tree base (downstream)
    q     ::FT = FT(0.0)    # mol s⁻¹ | flow rate in the xylem

    # pressure, k, and p_history profile
    k_element::Array{FT,1} =  ones(FT,10) .* FT(83.33)    # mol s⁻¹ MPa⁻¹ | a list of trunk k_max per element
    p_element::Array{FT,1} = zeros(FT,10)                 # MPa           | a list of trunk xylem pressure per element
    p_history::Array{FT,1} = zeros(FT,10)                 # MPa           | a list of trunk xylem pressure history per element
    t_element::Array{FT,1} =  ones(FT,10) .* K_25         # K             | a list of stem temperature for each element
    z_element::Array{FT,1} =  ones(FT,10) .* FT( 0.1 )    # m             | a list of trunk height per element
end




# function to create root list
function create_root_list(FT, n=5, z_lo=-1.0)
    # create a list of root layers
    root_list = [RootLayer{FT}() for i in 1:n]

    # iterate through the root list
    for i in 1:n
        # update the k_max for each layer, k_max = k_s * area / delta_h / n
        root_list[i].f_layer   = FT(1.0 / n)
        root_list[i].k_max     = FT(1.5625 / n)
        root_list[i].k_element = ones(FT,10) .* FT(15.625 / n)
        
        # update the root z, z_hi and z_lo for plotting purpose, z-element for computing gravitational pressure drop
        root_list[i].z_hi      = FT(z_lo * (i-1)/n)
        root_list[i].z_lo      = FT(z_lo * i/n)
        root_list[i].z_element = ones(FT,10) .* FT(z_lo * i/n * 0.1)
    end

    # return the root list
    return root_list
end




# struct for tree root
Base.@kwdef mutable struct Root{FT<:AbstractFloat,n}
    root_frac::Array{FT,1}            = ones(FT,5) .* FT(0.2)              # | a list of root fraction
    root_zs  ::Array{FT,1}            = FT.([-0.2,-0.4,-0.6,-0.8,-1.0])    # | a list of root z (lower value)
    root_list::Array{RootLayer{FT},1} = create_root_list(FT,n)             # | a list of root layer
end
