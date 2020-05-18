"""
    struct RootLayer{FT<:AbstractFloat}

A RootLayer type which contains root hydraulics information, including rhizosphere conductance.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct RootLayer{FT<:AbstractFloat}
    # root structure
    "Upper root depth `[m]`"
    z_hi   ::FT = FT( 0.0)
    "Lower root depth `[m]`"
    z_lo   ::FT = FT(-1.0)
    "Average horizental root length `[m]`"
    r_layer::FT = FT( 1.0)
    "Fraction of biomass and k_max"
    f_layer::FT = FT( 1.0)

    # hydraulic parameters
    "Weibull function (`k = k_max * exp( -(-p/B)^C )`) parameter B `[MPa]`"
    b     ::FT = FT( 2.0   )
    "Weibull function (`k = k_max * exp( -(-p/B)^C )`) parameter C"
    c     ::FT = FT( 5.0   )
    "Maximal hydraulic conductance, default at 2700.0 Kg h⁻¹ MPa⁻¹ m⁻² basal area `[mol s⁻¹ MPa⁻¹]`"
    k_max ::FT = FT(1.5625 )
    "Maximal hydraulic conductivity per cross section basal area per root depth `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_s   ::FT = FT(15.625 )
    "Maximal rhizosphere conductance, default at 9.72e12 Kg h⁻¹ MPa⁻¹ mm⁻² basal area `[mol s⁻¹ MPa⁻¹]`"
    k_rhiz::FT = FT( 1.5e10)

    # soil parameters
    "Soil water potential (upstream) `[MPa]`"
    p_ups    ::FT     = FT(  0.0   )
    "Soil texture parameter"
    soil_a   ::FT     = FT(602.0419)
    "Soil texture parameter, 1 - 1/n"
    soil_m   ::FT     = FT(  3.2432)
    "Soil texture parameter"
    soil_n   ::FT     = FT(  1.48  )
    "Maximal soil volumatric water content, vary with rock fraction"
    soil_msc ::FT     = FT(  0.4   )
    "Relative soil volumatric water content"
    soil_rwc ::FT     = FT(  1.0   )
    "Soil volumatric water content"
    soil_vwc ::FT     = FT(  0.4   )
    "Soil temperature `[K]`"
    t_soil   ::FT     = K_25
    "Soil texture"
    soil_type::String = "Sandy Clay Loam"

    # flows and pressures (need to be updated with time)
    "Water potential at the root-rhizosphere interface `[MPa]`"
    p_rhiz::FT = FT(0.0)
    "Xylem water pressure at the tree base (downstream) `[MPa]`"
    p_dos ::FT = FT(0.0)
    "Flow rate in the root layer `[mol s⁻¹]`"
    q     ::FT = FT(0.0)

    # pressure, k, and p_history profile
    "List of k_max per element `[mol s⁻¹ MPa⁻¹]`"
    k_element::Array{FT,1} =  ones(FT,10) .* FT(83.33)
    "List of xylem water pressure per element `[MPa]`"
    p_element::Array{FT,1} = zeros(FT,10)
    "List of xylem water pressure history (normalized to 298.15 K) per element `[MPa]`"
    p_history::Array{FT,1} = zeros(FT,10)
    "List of water temperature per element `[K]`"
    t_element::Array{FT,1} =  ones(FT,10) .* K_25
    "List of height per element to account for gravity `[m]`"
    z_element::Array{FT,1} =  ones(FT,10) .* FT( 0.1 )
end




"""
    create_root_list(FT, n, z_lo)

The function creates a root system with n root layers evenly in the soil. This function need to work with another function (`pending`) to quantify the root distribution.
- `FT` Floating type for the RootLayer struct
- `n` Total number of root layers in the root system
- `z_lo` Maximal root depth
"""
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




"""
    struct Root{FT, n}

A Root Struct which contains, by default, 5 even [`RootLayer`](@ref).

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct Root{FT<:AbstractFloat,n}
    "List of root fraction"
    root_frac::Array{FT,1}            = ones(FT,5) .* FT(0.2)
    "List of root z lower boundary"
    root_zs  ::Array{FT,1}            = FT.([-0.2,-0.4,-0.6,-0.8,-1.0])
    "List of [`RootLayer`](@ref)"
    root_list::Array{RootLayer{FT},1} = create_root_list(FT,n)
end
