# struct for leaf
Base.@kwdef mutable struct StructTreeLeaf
    # leaf structure
    angle_locat::Float32 = 0.0    # degree | direction from the tree trunk to the leaf, 0 for east, 270 for south, and 180 for west
    angle_incli::Float32 = 0.0    # degree | direction of leave, 0 for flat and 90 for vertical

    # leaf hydraulic parameters, per leaf area
    b     ::Float32 = 2.0     # MPa               | Weibull B
    c     ::Float32 = 5.0     #                   | Weibull C
    k_sla ::Float32 = 1.35    # mol s⁻¹ MPa⁻¹ m⁻² | maximal leaf hydraulic conductance per leaf area

    # flows and pressures (need to be updated with time)
    p_ups  ::Float32 = 0.0       # MPa          | xylem pressure at the leaf basa (upstream)
    p_i    ::Float32 = 0.0       # Pa           | leaf internal CO₂ partial pressure
    p_dos  ::Float32 = 0.0       # MPa          | xylem pressure of the leaf (downstream)

    # pressure, k, and p_history profile
    k_element::Array{Float32,1} = ones(10) * 13.5      # mol s⁻¹ MPa⁻¹ | a list of trunk k_max per element
    p_element::Array{Float32,1} = zeros(10)            # MPa           | a list of trunk xylem pressure per element
    p_history::Array{Float32,1} = zeros(10)            # MPa           | a list of trunk xylem pressure history per element
    t_element::Array{Float32,1} = ones(10) * 298.15    # K             | a list of stem temperature for each element
    z_element::Array{Float32,1} = ones(10) * 0.0       # m             | a list of trunk height per element
end




# struct of canopy layer
Base.@kwdef mutable struct StructTreeCanopyLayer
    # canopy layer and leaf structure
    f_layer::Float32 = 1.0      #    | fraction of LA of the layer, f_layer = 1.0 / n_layer
    f_view ::Float32 = 0.5      #    | view factor for sunlit leaves
    la     ::Float32 = 150.0    # m² | leaf area in the layer, la = la:ba * ba
    width  ::Float32 = 0.1      # m  | leaf width

    # canopy layer environemnt
    p_a  ::Float32 = 40.0        # Pa    | atmospheric CO₂ partial pressure
    p_atm::Float32 = 101325.0    # Pa    | atmospheric pressure
    p_O₂ ::Float32 = 21278.25    # Pa    | atmospheric O₂ partial pressure
    p_H₂O::Float32 = 1500.0      # Pa    | atmospheric H₂O partial pressure
    t_air::Float32 = 298.15      # K     | air temperature
    wind ::Float32 = 2.0         # m s⁻¹ | wind speed

    # leaf photosynthetic parameters
    g_ias_c::Float32 = 0.0      # unitless     | mesophyll conductance correction factor: multiplier
    g_ias_e::Float32 = 0.3      # unitless     | mesophyll conductance correction factor: exponent
    g_max  ::Float32 = 0.8      # mol m⁻² s⁻¹  | maximal gs at 25 °C
    Γ_star ::Float32 = 2.5      # Pa           | CO₂ compensation point with the absence of dark respiration
    gs_nssf::Float32 = 0.025    #              | non-steady state factor (use by multiplying the ∂A/∂E - ∂Θ/∂E)
    j_max  ::Float32 = 133.6    # μmol m⁻² s⁻¹ | maximal electron transport rate
    r_25   ::Float32 = 1.2      # μmol m⁻² s⁻¹ | leaf respiration rate
    v_max  ::Float32 = 80.0     # μmol m⁻² s⁻¹ | maximal carboxylation rate

    # leaf layers (e_list and q_list need to be updated with time)
    ag_list  ::Array{Float32,1}        = zeros(469)                           # μmol m⁻² s⁻¹ | gross a list per leaf area
    an_list  ::Array{Float32,1}        = zeros(469)                           # μmol m⁻² s⁻¹ | net a list per leaf area
    d_list   ::Array{Float32,1}        = zeros(469) .+ 0.015                  #              | leaf-to-air d list
    e_list   ::Array{Float32,1}        = zeros(469)                           # mol m⁻² s⁻¹  | flow rate list per leaf area
    ec_list  ::Array{Float32,1}        = zeros(469) .+ 6.3e-4                 # mol m⁻² s⁻¹  | e_crit list per leaf area
    gsc_list ::Array{Float32,1}        = zeros(469)                           # mol m⁻² s⁻¹  | gsc list per leaf area
    gsw_list ::Array{Float32,1}        = zeros(469)                           # mol m⁻² s⁻¹  | gsw list per leaf area
    la_list  ::Array{Float32,1}        = [ones(468)/468*75.0; 75.0]           # m²           | leaf area list
    leaf_list::Array{StructTreeLeaf,1} = [StructTreeLeaf() for i in 1:469]    #              | leaf struct list
    par_list ::Array{Float32,1}        = zeros(469)                           # μmol m⁻² s⁻¹ | PAR list
    pi_list  ::Array{Float32,1}        = zeros(469)                           # Pa           | leaf interanal CO₂ list
    q_list   ::Array{Float32,1}        = zeros(469)                           # mol s⁻¹      | flow rate list
    r_list   ::Array{Float32,1}        = zeros(469)                           # μmol m⁻² s⁻¹ | respiration list per leaf area
    t_list   ::Array{Float32,1}        = zeros(469) .+ 298.15                 # K            | leaf temperature list
end




# the struct for leaf
Base.@kwdef mutable struct StructTreeCanopy
    n_layer    ::Int                            = 20                                         # | number of canopy layers
    canopy_list::Array{StructTreeCanopyLayer,1} = [StructTreeCanopyLayer() for i in 1:20]    # | a list of leaf layers
end
