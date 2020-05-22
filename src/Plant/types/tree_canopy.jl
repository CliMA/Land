"""
    Leaf{FT<:AbstractFloat}

A Leaf type which contains leaf hydraulics information.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct Leaf{FT<:AbstractFloat}
    # leaf angles, may delete this in the future
    "direction from the tree trunk to the leaf, 0 for east, 270 for south, and 180 for west `[°]`"
    angle_locat::FT = FT(0.0)
    "direction of leave, 0 for flat and 90 for vertical `[°]`"
    angle_incli::FT = FT(0.0)

    # leaf hydraulic parameters
    "Weibull function (`k = k_max * exp( -(-p/B)^C )`) parameter B `[MPa]`"
    b    ::FT = FT(2.0 )
    "Weibull function (`k = k_max * exp( -(-p/B)^C )`) parameter C"
    c    ::FT = FT(5.0 )
    "Maximal leaf hydraulic conductance per leaf area `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_sla::FT = FT(1.35)

    # flows and pressures (need to be updated with time)
    "Leaf xylem water pressure (different from leaf water potential) at the leaf base (upstream) `[MPa]`"
    p_ups::FT = FT(0.0)
    "Leaf xylem water pressure at the downstream end of leaf xylem `[MPa]`"
    p_dos::FT = FT(0.0)

    # pressure, k, and p_history profile
    "List of leaf k_max per element (mol s⁻¹ MPa⁻¹)"
    k_element::Array{FT,1} =  ones(FT,10) * FT( 13.5 )
    "List of xylem water pressure per element `[MPa]`"
    p_element::Array{FT,1} = zeros(FT,10)
    "List of xylem water pressure history (normalized to 298.15 K) per element `[MPa]`"
    p_history::Array{FT,1} = zeros(FT,10)
    "List of xylem water temperature per element `[K]`"
    t_element::Array{FT,1} =  ones(FT,10) * FT(298.15)
    "List of leaf element height change to account for gravity, default = 0 `[m]`"
    z_element::Array{FT,1} =  ones(FT,10) * FT(  0.0 )
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
    f_layer::FT = FT(  1.0)
    "Fraction of the sunlit leaves"
    f_view ::FT = FT(  0.5)
    "Leaf area in the layer `[m²]`"
    la     ::FT = FT(150.0)
    "Leaf width `[m]`"
    width  ::FT = FT(  0.1)

    # canopy layer environemnt
    "Atmospheric CO₂ partial pressure in the layer `[Pa]`"
    p_a  ::FT = FT(    40.0 )
    "Atmospheric pressure `[Pa]`"
    p_atm::FT = FT(101325.0 )
    "Atmospheric O₂ partial presssure `[Pa]`"
    p_O₂ ::FT = FT( 21278.25)
    "Atmospheric H₂O partial pressure `[Pa]`"
    p_H₂O::FT = FT(  1500.0 )
    "Air temperature `[K]`"
    t_air::FT = FT(   298.15)
    "Wind speed `[m s⁻¹]`"
    wind ::FT = FT(     2.0 )

    # leaf photosynthetic parameters
    "Mesophyll conductance correction factor: multiplier"
    g_ias_c  ::FT = FT(  0.0  )
    "Mesophyll conductance correction factor: exponent"
    g_ias_e  ::FT = FT(  0.3  )
    "Maximal leaf diffusive conductance for H₂O at 298.15 K `[mol m⁻² s⁻¹]`"
    g_max    ::FT = FT(  0.8  )
    "CO₂ compensation point with the absence of dark respiration `[Pa]`"
    Γ_star   ::FT = FT(  2.5  )
    "Non-steady state factor (use it by multiplying the ∂A/∂E - ∂Θ/∂E)"
    gs_nssf  ::FT = FT(  0.025)
    "Maximal electron transport rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    j_max    ::FT = FT(133.6  )
    "Leaf respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_25     ::FT = FT(  1.2  )
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_max    ::FT = FT( 80.0  )
    "Maximal PEP carboxylation rate at 298.15 K (only for C4 plants) `[μmol m⁻² s⁻¹]`"
    p_max    ::FT = FT(120.0  )
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
    "List of leaf-to-air vapor pressure deficit `[unitless]`"
    d_list   ::Array{FT,1} = zeros(FT,n_leaf) .+ FT(0.015)
    "List of flow rate per leaf area `[mol m⁻² s⁻¹]`"
    e_list   ::Array{FT,1} = zeros(FT,n_leaf)
    "List of critical flow rate (when leaf xylem pressure induces desiccation) per leaf area `[mol m⁻² s⁻¹]`"
    ec_list  ::Array{FT,1} = zeros(FT,n_leaf) .+ FT(6.3e-4)
    "List of effective leaf diffusive conductance for CO₂ `[mol m⁻² s⁻¹]`"
    gsc_list ::Array{FT,1} = zeros(FT,n_leaf)
    "List of leaf diffusive conductance for H₂O `[mol m⁻² s⁻¹]`"
    gsw_list ::Array{FT,1} = zeros(FT,n_leaf)
    "List of leaf area per leaf `[m²]`"
    la_list  ::Array{FT,1} = FT.([ones(n_leaf-1)/(n_leaf-1)*75.0; 75.0])
    "List of PAR for each leaf"
    par_list ::Array{FT,1} = zeros(FT,n_leaf)
    "List of leaf internal CO₂ partial pressure `[Pa]`"
    pi_list  ::Array{FT,1} = zeros(FT,n_leaf)
    "List of flow rate `[mol s⁻¹]`"
    q_list   ::Array{FT,1} = zeros(FT,n_leaf)
    "List of leaf respiration rate `[μmol m⁻² s⁻¹]`"
    r_list   ::Array{FT,1} = zeros(FT,n_leaf)
    "List of leaf temperature `[K]`"
    t_list   ::Array{FT,1} = zeros(FT,n_leaf) .+ FT(298.15)
end

"""
    Canopy{FT<:AbstractFloat,n, n_leaf}

A Canopy type which contains `n` [`CanopyLayer`](@ref).

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct Canopy{FT<:AbstractFloat, n, n_leaf}
    "Canopy layer number"
    n_layer::Int = n
    "List of [`CanopyLayer`](@ref)"
    canopy_list::Array{CanopyLayer{FT},1} = [CanopyLayer{FT,n_leaf}() for i in 1:n]
end
