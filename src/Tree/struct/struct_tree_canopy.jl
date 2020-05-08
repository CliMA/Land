# struct for leaf
Base.@kwdef mutable struct struct_tree_leaf
    # leaf structure
    angle_locat::Float32 = 0.0    # degree | direction from the tree trunk to the leaf, 0 for east, 270 for south, and 180 for west
    angle_incli::Float32 = 0.0    # degree | direction of leave, 0 for flat and 90 for vertical
    width      ::Float32 = 0.1    # m      | leaf width

    # leaf hydraulic parameters, per leaf area
    b     ::Float32 = 2.0     # MPa                  | Weibull B
    c     ::Float32 = 5.0     #                      | Weibull C
    k_sla ::Float32 = 1.35    # mol s^-1 MPa^-1 m^-2 | maximal leaf hydraulic conductance per leaf area

    # leaf photosynthetic parameters
    g_max::Float32 = 0.8      # mol m^-2 s^-1  | maximal gs at 25 degree C
    j_max::Float32 = 133.6    # umol m^-2 s^-1 | maximal electron transport rate
    par  ::Float32 = 0.0      # umol m^-2 s^-1 | photosynthetic active radiation
    v_max::Float32 = 80.0     # umol m^-2 s^-1 | maximal carboxylation rate

    # flows and pressures
    a_gross::Float32 = 0.0       # umol m^-2 s^-1 | gross photosynthetic rate
    a_net  ::Float32 = 0.0       # umol m^-2 s^-1 | net photosynthetic rate
    e      ::Float32 = 0.0       # mol m^-2 s^-1  | flow rate in the xylem
    gsc    ::Float32 = 0.0       # mol m^-2 s^-1  | stomatal conductance for CO2
    gsw    ::Float32 = 0.0       # mol m^-2 s^-1  | stomatal conductance for H2O
    gs_nssf::Float32 = 1e2       #                | non-steady state factor (use by multiplying the d_optimizer / d_g)
    p_base ::Float32 = 0.0       # MPa            | xylem pressure at the leaf basa
    p_i    ::Float32 = 0.0       # Pa             | leaf internal CO2
    p_leaf ::Float32 = 0.0       # MPa            | xylem pressure of the leaf
    r      ::Float32 = 0.0       # umol m^-2 s^-1 | respiration rate
    t_leaf ::Float32 = 198.15    # K              | leaf temperature

    # pressure, k, and p_history profile
    k_element::Array{Float32,1} = ones(10) * 13.5      # mol s^-1 MPa^-1 | a list of trunk k_max per element
    p_element::Array{Float32,1} = zeros(10)            # MPa             | a list of trunk xylem pressure per element
    p_history::Array{Float32,1} = zeros(10)            # MPa             | a list of trunk xylem pressure history per element
    t_element::Array{Float32,1} = ones(10) * 298.15    # K               | a list of stem temperature for each element
    z_element::Array{Float32,1} = ones(10) * 0.0       # m               | a list of trunk height per element
end




# struct of canopy layer
Base.@kwdef mutable struct struct_tree_canopy_layer
    # canopy layer structure
    f_layer::Float32 = 1.0      #     | fraction of LA of the layer, f_layer = 1.0 / n_layer
    f_view ::Float32 = 0.5      #     | view factor for sunlit leaves
    la     ::Float32 = 150.0    # m^2 | leaf area in the layer, la = la:ba * ba

    # canopy layer environemnt
    p_a  ::Float32 = 40.0        # Pa     | atmospheric CO2 partial pressure
    p_atm::Float32 = 101325.0    # Pa     | atmospheric pressure
    p_o2 ::Float32 = 21278.25    # Pa     | atmospheric O2 partial pressure
    p_h2o::Float32 = 1500.0      # Pa     | atmospheric H2O partial pressure
    t_air::Float32 = 298.15      # K      | air temperature
    wind ::Float32 = 2.0         # m s^-1 | wind speed

    # leaf layers
    e_list   ::Array{Float32,1}          = zeros(469)                             # mol m^-2 s^1   | flow rate list per leaf area
    la_list  ::Array{Float32,1}          = [ones(468)/468*75.0; 75.0]             # m^2            | leaf area list
    leaf_list::Array{struct_tree_leaf,1} = [struct_tree_leaf() for i in 1:469]    #                | leaf struct list
    par_list ::Array{Float32,1}          = zeros(469)                             # umol m^-2 s^-1 | PAR list
    q_list   ::Array{Float32,1}          = zeros(469)                             # mol s^1        | flow rate list
end




# the struct for leaf
Base.@kwdef mutable struct struct_tree_canopy
    canopy_list ::Array{struct_tree_canopy_layer,1} = [struct_tree_canopy_layer() for i in 1:20]    # | a list of leaf layers
end