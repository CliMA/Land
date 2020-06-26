#=
###############################################################################
#
# Root type and function to create root system
# This struct and function passed the FT test
# This struct and function are documented in the Plant page
#
###############################################################################
"""
    RootLayer{FT<:AbstractFloat}

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
    "Maximal rhizosphere conductance, default at 9.72e12 Kg h⁻¹ MPa⁻¹ m⁻² basal area `[mol s⁻¹ MPa⁻¹]`"
    k_rhiz::FT = FT( 1.5e10)

    # soil parameters
    "Soil water potential (upstream) `[MPa]`"
    p_ups    ::FT     = FT(  0.0   )
    "Soil texture parameter"
    soil_a   ::FT     = FT(602.0419)
    "Soil texture parameter, 1 - 1/n"
    soil_m   ::FT     = FT( 0.32432)
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

    return root_list
end








###############################################################################
#
# Stem type and function to create branch
# These structs passed the FT test
# This struct and function are documented in the Plant page
#
###############################################################################
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








###############################################################################
#
# Canopy and leaf types and function to create branch
# These structs passed the FT test
# These structs are documented in the Plant page
#
###############################################################################
"""
    Leaf{FT<:AbstractFloat}

A Leaf type which contains leaf hydraulics information.

Migrated to Hydraulics module

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct Leaf{FT<:AbstractFloat}
    # leaf angles, may delete this in the future
    "direction from the tree trunk to the leaf, 0 for east, 270 for south, and 180 for west `[°]`"
    angle_locat::FT = FT(0.0)
    "direction of leave, 0 for flat and 90 for vertical `[°]`"
    angle_incli::FT = FT(0.0)
end




"""
    CanopyLayer{FT<:AbstractFloat, n_leaf}

A CanopyLayer type, which constains environemntal conditions and leaf-level fluxes for `n_leaf` [`Leaf`](@ref). The `n_leaf` is the sum of `n_Ari * n_Incli` sunlit leaves and 1 shaded leaves.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct CanopyLayer{FT<:AbstractFloat, n_leaf}
    # canopy layer and leaf structure
    "Leaf area fraction in the canopy"
    f_layer::FT = FT(1.0)
    "Fraction of the sunlit leaves"
    f_view ::FT = FT(0.5)
    "Leaf area in the layer `[m²]`"
    la     ::FT = FT(7.5)
    "Leaf width `[m]`"
    width  ::FT = FT(0.1)

    # canopy layer environemnt
    "Atmospheric CO₂ partial pressure in the layer `[Pa]`"
    p_a  ::FT = FT(40.0)
    "Atmospheric pressure `[Pa]`"
    p_atm::FT = FT(101325.0)
    "Atmospheric O₂ partial presssure `[Pa]`"
    p_O₂ ::FT = FT(21278.25)
    "Atmospheric H₂O partial pressure `[Pa]`"
    p_H₂O::FT = FT(1500.0)
    "Air temperature `[K]`"
    t_air::FT = FT(298.15)
    "Wind speed `[m s⁻¹]`"
    wind ::FT = FT(2.0)

    # leaf photosynthetic parameters
    "Mesophyll conductance correction factor: multiplier"
    g_ias_c  ::FT = FT(0.0)
    "Mesophyll conductance correction factor: exponent"
    g_ias_e  ::FT = FT(0.3)
    "Maximal leaf diffusive conductance for H₂O at 298.15 K `[mol m⁻² s⁻¹]`"
    g_max    ::FT = FT(0.8)
    "Minimal leaf diffusive conductance for H₂O at 298.15 K `[mol m⁻² s⁻¹]`"
    g_min    ::FT = FT(0.0016)
    "Non-steady state factor for optimization approach (use it by multiplying the ∂A/∂E - ∂Θ/∂E)"
    gs_nssf  ::FT = FT(0.025)
    "Non-steady state factor for empirical apporach (used by multiplying the gsw_mod - gsw_curr)"
    gs_empi  ::FT = FT(0.003)
    "Maximal electron transport rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    j_max    ::FT = FT(133.6)
    "Leaf respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_25     ::FT = FT(1.2)
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_max    ::FT = FT(80.0)
    "Maximal PEP carboxylation rate at 298.15 K (only for C4 plants) `[μmol m⁻² s⁻¹]`"
    p_max    ::FT = FT(120.0)
    "Curvature ratio for J"
    curvature::FT = FT(0.9)
    "Quntum yield for electron, using aPAR"
    qy       ::FT = FT(0.4081632653061224)

    # leaf layers (e_list and q_list need to be updated with time)
    "List of [`Leaf`](@ref)"
    leaf_list::Array{Leaf{FT},1} = [Leaf{FT}() for i in 1:n_leaf]
    "List of gross photosynthetic rate for n_leaf leaves `[μmol m⁻² s⁻¹]`"
    ag_list  ::Array{FT,1} = zeros(FT,n_leaf)
    "List of net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    an_list  ::Array{FT,1} = zeros(FT,n_leaf)
    "List if maximal A at the given scenario `[μmol m⁻² s⁻¹]`"
    am_list  ::Array{FT,1} = zeros(FT,n_leaf)
    "List of leaf-to-air vapor pressure deficit `[unitless]`"
    d_list   ::Array{FT,1} = zeros(FT,n_leaf) .+ FT(0.015)
    "List of flow rate per leaf area `[mol m⁻² s⁻¹]`"
    e_list   ::Array{FT,1} = zeros(FT,n_leaf)
    "List of critical flow rate (when leaf xylem pressure induces desiccation) per leaf area `[mol m⁻² s⁻¹]`"
    ec_list  ::Array{FT,1} = zeros(FT,n_leaf) .+ FT(0.012467781)
    "List of effective leaf diffusive conductance for CO₂ `[mol m⁻² s⁻¹]`"
    glc_list ::Array{FT,1} = zeros(FT,n_leaf) .+ FT(0.001)
    "List of empirical leaf diffusive conductance for H₂O `[mol m⁻² s⁻¹]`"
    glw_empi ::Array{FT,1} = zeros(FT,n_leaf) .+ FT(0.0016)
    "List of leaf diffusive conductance for H₂O `[mol m⁻² s⁻¹]`"
    glw_list ::Array{FT,1} = zeros(FT,n_leaf) .+ FT(0.0016)
    "List of maximal hydraulic conductance [mol m⁻² s⁻¹]"
    km_list  ::Array{FT,1} =  ones(FT,n_leaf)
    "List of leaf area per leaf `[m²]`"
    la_list  ::Array{FT,1} = FT.([ones(n_leaf-1)/(n_leaf-1)*3.75; 3.75])
    "List of PAR for each leaf `[μmol m⁻² s⁻¹]`"
    par_list ::Array{FT,1} = zeros(FT,n_leaf) .+ 100
    "List of leaf internal CO₂ partial pressure `[Pa]`"
    pi_list  ::Array{FT,1} = zeros(FT,n_leaf)
    "List of flow rate `[mol s⁻¹]`"
    q_list   ::Array{FT,1} = zeros(FT,n_leaf)
    "List of leaf respiration rate `[μmol m⁻² s⁻¹]`"
    r_list   ::Array{FT,1} = zeros(FT,n_leaf)
    "List of leaf temperature `[K]`"
    t_list   ::Array{FT,1} = zeros(FT,n_leaf) .+ FT(298.15)
end








###############################################################################
#
# Canopy type and function to create branch
# This struct passed the FT test
# This struct is documented in the Plant page
#
###############################################################################
"""
    Tree{FT<:AbstractFloat}

A Tree type which includes
- `root_no` number of root layers
- a trunk
- `canopy_no` branches
- `leaf_no` canopy layers

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct Tree{FT<:AbstractFloat, root_no, canopy_no, leaf_no}
    # plant information
    "Age `[year]`"
    age::Int = 10
    "Basal area `[m²]`"
    ba ::FT  = FT( 0.1)
    "Ground area `[m²]`"
    ga ::FT  = FT(50.0)
    "Tree height `[m]`"
    h  ::FT  = FT( 8.0)
    "Photosynthesis type C3/C4/CAM"
    photo_type::String = "C3"
    "Photosynthesis model parameter set, type undefined"
    photo_para_set::AbstractPhotoModelParaSet = C3Bernacchi(FT)
    "Optimization model option, type undefined"
    stomata_scheme::AbstractStomatalModel = OSMWang()

    # TODO change the code accordingly, particularly the optimization models
    "Use steady-state or non-steady-state stomatal behavior"
    steady_state::Bool = false

    # tree formation of root
    "Number of root layers"
    n_root   ::Int                    = root_no
    "List of root fraction"
    root_frac::Array{FT,1}            = ones(FT,5) .* FT(0.2)
    "List of root z lower boundary"
    root_zs  ::Array{FT,1}            = FT.([-0.2,-0.4,-0.6,-0.8,-1.0])
    "List of [`RootLayer`](@ref)"
    root_list::Array{RootLayer{FT},1} = create_root_list(FT,n_root)

    # tree information of trunk and branch
    "Number of canopy layers"
    n_canopy   ::Int               = canopy_no
    "Trunk using type [`Stem`](@ref)"
    trunk      ::Stem              = Stem{FT}()
    "List of [`Stem`](@ref)"
    branch_list::Array{Stem{FT},1} = create_branch_list(FT,n_canopy)

    # tree information of canopy
    "Number of leaves in each canopy layer"
    n_leaf     ::Int                      = leaf_no
    "List of [`CanopyLayer`](@ref)"
    canopy_list::Array{CanopyLayer{FT},1} = [CanopyLayer{FT,n_leaf}() for i in 1:n_canopy]
end
=#