# struct for leaf, may need to merge with the leaf struct in the RT module
Base.@kwdef mutable struct Leaf{FT<:AbstractFloat}
    # leaf structure
    angle_locat::FT = FT(0.0)    # degree | direction from the tree trunk to the leaf, 0 for east, 270 for south, and 180 for west
    angle_incli::FT = FT(0.0)    # degree | direction of leave, 0 for flat and 90 for vertical

    # leaf hydraulic parameters, per leaf area
    b    ::FT = FT(2.0 )    # MPa               | Weibull B
    c    ::FT = FT(5.0 )    #                   | Weibull C
    k_sla::FT = FT(1.35)    # mol s⁻¹ MPa⁻¹ m⁻² | maximal leaf hydraulic conductance per leaf area

    # flows and pressures (need to be updated with time)
    p_ups::FT = FT(0.0)    # MPa          | xylem pressure at the leaf basa (upstream)
    p_i  ::FT = FT(0.0)    # Pa           | leaf internal CO₂ partial pressure
    p_dos::FT = FT(0.0)    # MPa          | xylem pressure of the leaf (downstream)

    # pressure, k, and p_history profile
    k_element::Array{FT,1} =  ones(FT,10) * FT( 13.5 )    # mol s⁻¹ MPa⁻¹ | a list of trunk k_max per element
    p_element::Array{FT,1} = zeros(FT,10)                 # MPa           | a list of trunk xylem pressure per element
    p_history::Array{FT,1} = zeros(FT,10)                 # MPa           | a list of trunk xylem pressure history per element
    t_element::Array{FT,1} =  ones(FT,10) * FT(298.15)    # K             | a list of stem temperature for each element
    z_element::Array{FT,1} =  ones(FT,10) * FT(  0.0 )    # m             | a list of trunk height per element
end




# struct of canopy layer
Base.@kwdef mutable struct CanopyLayer{FT<:AbstractFloat, n_total}
    # canopy layer and leaf structure
    f_layer::FT = FT(  1.0)    #    | fraction of LA of the layer, f_layer = 1.0 / n_layer
    f_view ::FT = FT(  0.5)    #    | view factor for sunlit leaves
    la     ::FT = FT(150.0)    # m² | leaf area in the layer, la = la:ba * ba
    width  ::FT = FT(  0.1)    # m  | leaf width

    # canopy layer environemnt
    p_a  ::FT = FT(    40.0 )    # Pa    | atmospheric CO₂ partial pressure
    p_atm::FT = FT(101325.0 )    # Pa    | atmospheric pressure
    p_O₂ ::FT = FT( 21278.25)    # Pa    | atmospheric O₂ partial pressure
    p_H₂O::FT = FT(  1500.0 )    # Pa    | atmospheric H₂O partial pressure
    t_air::FT = FT(   298.15)    # K     | air temperature
    wind ::FT = FT(     2.0 )    # m s⁻¹ | wind speed

    # leaf photosynthetic parameters
    g_ias_c::FT = FT(  0.0  )    # unitless     | mesophyll conductance correction factor: multiplier
    g_ias_e::FT = FT(  0.3  )    # unitless     | mesophyll conductance correction factor: exponent
    g_max  ::FT = FT(  0.8  )    # mol m⁻² s⁻¹  | maximal gs at 25 °C
    Γ_star ::FT = FT(  2.5  )    # Pa           | CO₂ compensation point with the absence of dark respiration
    gs_nssf::FT = FT(  0.025)    #              | non-steady state factor (use by multiplying the ∂A/∂E - ∂Θ/∂E)
    j_max  ::FT = FT(133.6  )    # μmol m⁻² s⁻¹ | maximal electron transport rate
    r_25   ::FT = FT(  1.2  )    # μmol m⁻² s⁻¹ | leaf respiration rate
    v_max  ::FT = FT( 80.0  )    # μmol m⁻² s⁻¹ | maximal carboxylation rate

    # leaf layers (e_list and q_list need to be updated with time)
    ag_list  ::Array{FT,1}       = zeros(FT,n_total)                                # μmol m⁻² s⁻¹ | gross a list per leaf area
    an_list  ::Array{FT,1}       = zeros(FT,n_total)                                # μmol m⁻² s⁻¹ | net a list per leaf area
    d_list   ::Array{FT,1}       = zeros(FT,n_total) .+ FT(0.015)                   #              | leaf-to-air d list
    e_list   ::Array{FT,1}       = zeros(FT,n_total)                                # mol m⁻² s⁻¹  | flow rate list per leaf area
    ec_list  ::Array{FT,1}       = zeros(FT,n_total) .+ FT(6.3e-4)                  # mol m⁻² s⁻¹  | e_crit list per leaf area
    gsc_list ::Array{FT,1}       = zeros(FT,n_total)                                # mol m⁻² s⁻¹  | gsc list per leaf area
    gsw_list ::Array{FT,1}       = zeros(FT,n_total)                                # mol m⁻² s⁻¹  | gsw list per leaf area
    la_list  ::Array{FT,1}       = FT.([ones(n_total-1)/(n_total-1)*75.0; 75.0])    # m²           | leaf area list
    leaf_list::Array{Leaf{FT},1} = [Leaf{FT}() for i in 1:n_total]                  #              | leaf struct list
    par_list ::Array{FT,1}       = zeros(FT,n_total)                                # μmol m⁻² s⁻¹ | PAR list
    pi_list  ::Array{FT,1}       = zeros(FT,n_total)                                # Pa           | leaf interanal CO₂ list
    q_list   ::Array{FT,1}       = zeros(FT,n_total)                                # mol s⁻¹      | flow rate list
    r_list   ::Array{FT,1}       = zeros(FT,n_total)                                # μmol m⁻² s⁻¹ | respiration list per leaf area
    t_list   ::Array{FT,1}       = zeros(FT,n_total) .+ FT(298.15 )                 # K            | leaf temperature list
end




# the struct for canopy, may need to merge with the canopy struct in the RT module
Base.@kwdef mutable struct Canopy{FT<:AbstractFloat,n, n_total}
    n_layer    ::Int                      = n                                           # | number of canopy layers
    canopy_list::Array{CanopyLayer{FT},1} = [CanopyLayer{FT,n_total}() for i in 1:n]    # | a list of leaf layers
end
