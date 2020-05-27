"""
    Stem{FT<:AbstractFloat}

A Stem type which contains hydraulic information, including drought history profile.
The Stem type is used for creating a trunk and a list of branches.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct Stem{FT<:AbstractFloat}
    # trunk structure
    "Horizental radius of the layer `[m]`"
    r_layer::FT = FT(1.0)
    "Upper height from the ground (trunk-branch joint for trunk, stem-leaf joint for branch) `[m]`"
    z_hi   ::FT = FT(3.0)
    "Lower height from the ground (0 for trunk, trunk-branch joint for branch) `[m]`"
    z_lo   ::FT = FT(0.0)

    # hydraulic parameter
    "Weibull function (`k = k_max * exp( -(-p/B)^C )`) parameter B `[MPa]`"
    b    ::FT = FT(  2.0  )
    "Weibull function (`k = k_max * exp( -(-p/B)^C )`) parameter C"
    c    ::FT = FT(  5.0  )
    "Maximal hydraulic conductance, = 5400.0 * 8/3 Kg h⁻¹ MPa⁻¹ m⁻² basal area `[mol s⁻¹ MPa⁻¹]`"
    k_max::FT = FT(  8.333)
    "Maximal hydraulic conductivity per cross section area per tree height `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_s  ::FT = FT(250.0  )

    # flows and pressures (need to be updated with time)
    "Xylem water pressure at the tree basa (upstream) `[MPa]`"
    p_ups::FT = FT(0.0)
    "Xylem water pressure at the trunk-stem or stem-leaf joint (downstream) `[MPa]`"
    p_dos::FT = FT(0.0)
    "Flow rate in the xylem `[mol s⁻¹]`"
    q    ::FT = FT(0.0)

    # pressure, k, and p_history profile
    "List of stem k_max per element `[mol s⁻¹ MPa⁻¹]`"
    k_element::Array{FT,1} =  ones(FT,10) .* FT( 83.33)
    "List of stem xylem pressure per element `[MPa]`"
    p_element::Array{FT,1} = zeros(FT,10)
    "List of stem xylem pressure history (normalized to 298.15 K) per element `[MPa]`"
    p_history::Array{FT,1} = zeros(FT,10)
    "List of stem temperature for each element `[K]`"
    t_element::Array{FT,1} =  ones(FT,10) .* FT(298.15)
    "List of stem height per element to account for gravity `[m]`"
    z_element::Array{FT,1} =  ones(FT,10) .* FT(  0.3 )
end




"""
    create_branch_list(FT, n, z_lo, z_hi)

This function creates a list of stem structs as branch system for the tree, given
- `FT` Floating type for the Stem struct
- `n` Total number of canopy layers, one branch for a canopy layer
- `z_lo` Minimal canopy distance from the ground
- `z_hi` Maximal canopy distance from the ground
Each branch (Stem) is responsible for a canopy layer. The branch index from 1:n is from lower to higher canopy.
"""
function create_branch_list(FT, n=20, z_lo=3.0, z_hi=8.0)
    # create a list of branches following the default of stem
    branch_list = [Stem{FT}() for i in 1:n]

    # iterate through the branch for each canopy layers
    for i in 1:n
        # update the k_max for each layer, k_max = k_s * area / delta_h / n
        branch_list[i].k_max     = FT(5.0 / n)
        branch_list[i].k_element = ones(FT,10) * FT(50.0 / n)
        
        # update the z for each layer, the z_lo and z_hi are meant for plotting, the z_element is meant for computing pressure drop
        branch_list[i].z_lo      = FT(z_lo + (z_hi-z_lo) * (i-1)/n)
        branch_list[i].z_hi      = FT(z_lo + (z_hi-z_lo) * i/n)
        branch_list[i].z_element = ones(FT,10) .* FT((z_hi-z_lo) * i/n * 0.1)
    end

    return branch_list
end
