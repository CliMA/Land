#=
###############################################################################
#
# All plant parameters in one struct
#
###############################################################################
Base.@kwdef mutable struct PlantParaSet{FT<:AbstractFloat}
    # Basic information
    "Basal area `[m²]`"
    ba::FT  = FT(0.1)
    "Ground area `[m²]`"
    ga::FT  = FT(50)
    "Leaf area index"
    lai::FT = FT(3)
    "Canopy layers"
    n_canopy::Int = 20
    "Leaf numbers in a canopy layer"
    n_leaf  ::Int = 325
    "Root layers"
    n_root  ::Int = 5

    # Rhizosphere hydraulics (5 layers by default)
    "Soil water potential"
    p_soil ::Array{FT,1} = zeros(FT,5)
    "Rhizosphere hydraulic conductance `[mol s⁻¹ MPa⁻¹]`"
    rhiz_k ::Array{FT,2} =  ones(FT,(5,10)) .* 5e10
    "Rhizosphere maximal hydraulic conductance `[mol s⁻¹ MPa⁻¹]`"
    rhiz_km::Array{FT,2} =  ones(FT,(5,10)) .* 5e10
    "Soil texture parameter"
    soil_a::FT = FT(602.0419)
    "Soil texture parameter, 1 - 1/n"
    soil_m::FT = FT(0.324324)
    "Soil texture parameter"
    soil_n::FT = FT(1.48)

    # Root hydraulics (5 layers by default)
    "Root Weibull B `[MPa]`"
    root_b ::Array{FT,1} =  ones(FT,5) .* 2
    "Root Weibull C"
    root_c ::Array{FT,1} =  ones(FT,5) .* 6
    "Root hydraulic conductance `[mol s⁻¹ MPa⁻¹]`"
    root_k ::Array{FT,2} =  ones(FT,(5,10)) .* 4
    "Root maximal hydraulic conductance `[mol s⁻¹ MPa⁻¹]`"
    root_km::Array{FT,2} =  ones(FT,(5,10)) .* 4
    "Root xylem pressure `[MPa]`"
    root_p ::Array{FT,2} = zeros(FT,(5,10))
    "Root element gravity pressure drop `[MPa]`"
    root_pg::Array{FT,2} = zeros(FT,(5,10))
    "Root xylem pressure history (normalized to 298.15 K) `[MPa]`"
    root_ph::Array{FT,2} = zeros(FT,(5,10))
    "Root xylem water relative surface tension"
    root_st::Array{FT,1} =  ones(FT,5)
    "Root xylem water temperature in each layer `[K]`"
    root_t ::Array{FT,1} = zeros(FT,5) .+ K_25
    "Root xylem water relative viscosity"
    root_vs::Array{FT,1} =  ones(FT,5)

    # Trunk hydraulics
    "Root Weibull B `[MPa]`"
    trunk_b ::FT = 2
    "Root Weibull C"
    trunk_c ::FT = 6
    "Trunk hydraulic conductance `[mol s⁻¹ MPa⁻¹]`"
    trunk_k ::Array{FT,1} =  ones(FT,10) .* 100
    "Trunk maximal hydraulic conductance `[mol s⁻¹ MPa⁻¹]`"
    trunk_km::Array{FT,1} =  ones(FT,10) .* 100
    "Trunk xylem pressure `[MPa]`"
    trunk_p ::Array{FT,1} = zeros(FT,10)
    "Trunk element gravity pressure drop `[MPa]`"
    trunk_pg::Array{FT,1} = zeros(FT,10)
    "Trunk xylem pressure history (normalized to 298.15 K) `[MPa]`"
    trunk_ph::Array{FT,1} = zeros(FT,10)
    "Trunk xylem water relative surface tension"
    trunk_st::FT = FT(1)
    "Trunk xylem water temperature in each layer `[K]`"
    trunk_t ::FT = FT(K_25)
    "Trunk xylem water relative viscosity"
    trunk_vs::FT = FT(1)

    # Stem hydraulics (20 layers by default)
    "Stem Weibull B `[MPa]`"
    stem_b ::FT = 2
    "Stem Weibull C"
    stem_c ::FT = 6
    "Stem hydraulic conductance `[mol s⁻¹ MPa⁻¹]`"
    stem_k ::Array{FT,2} =  ones(FT,(20,10)) .* 6
    "Stem maximal hydraulic conductance `[mol s⁻¹ MPa⁻¹]`"
    stem_km::Array{FT,2} =  ones(FT,(20,10)) .* 6
    "Stem xylem pressure `[MPa]`"
    stem_p ::Array{FT,2} = zeros(FT,(20,10))
    "Stem element gravity pressure drop `[MPa]`"
    stem_pg::Array{FT,2} = zeros(FT,(20,10))
    "Stem xylem pressure history (normalized to 298.15 K) `[MPa]`"
    stem_ph::Array{FT,2} = zeros(FT,(20,10))
    "Stem xylem water relative surface tension"
    stem_st::FT = FT(1)
    "Stem xylem water temperature in each layer `[K]`"
    stem_t ::FT = FT(K_25)
    "Stem xylem water relative viscosity"
    stem_vs::FT = FT(1)

    # Leaf hydraulics (20 layers of canopy, and 325 leaves per layer)
    "Leaf Weibull B `[MPa]`"
    leaf_b ::FT = 2
    "Leaf Weibull C"
    leaf_c ::FT = 6
    "Leaf hydraulic conductance `[mol s⁻¹ MPa⁻¹]`"
    leaf_k ::Array{FT,3} =  ones(FT,(325,10,20)) .* 1.5
    "Leaf maximal hydraulic conductance `[mol s⁻¹ MPa⁻¹]`"
    leaf_km::Array{FT,3} =  ones(FT,(325,10,20)) .* 1.5
    "Leaf area for each element `[m²]`"
    leaf_la::Array{FT,2} = zeros(FT,325,20)
    "Leaf xylem pressure `[MPa]`"
    leaf_p ::Array{FT,3} = zeros(FT,(325,10,20))
    "Leaf element gravity pressure drop `[MPa]`"
    leaf_pg::Array{FT,3} = zeros(FT,(325,10,20))
    "Leaf xylem pressure history (normalized to 298.15 K) `[MPa]`"
    leaf_ph::Array{FT,3} = zeros(FT,(325,10,20))
    "Leaf xylem water relative surface tension"
    leaf_st::Array{FT,2} =  ones(FT,325,20)
    "Leaf xylem water temperature in each layer `[K]`"
    leaf_t ::Array{FT,2} = zeros(FT,325,20) .+ K_25
    "Leaf xylem water relative viscosity"
    leaf_vs::Array{FT,2} =  ones(FT,325,20)

    # Leaf gas exchange --- stomata
    "Transpiration rate per leaf area `[mol m⁻² s⁻¹]`"
    leaf_e   ::Array{FT,2} = zeros(FT,(325,20))
    "Critical transpiration rate per leaf area `[mol m⁻² s⁻¹]`"
    leaf_ec  ::Array{FT,2} = zeros(FT,(325,20))
    "Boundary layer conductance to H₂O `[mol m⁻² s⁻¹]`"
    leaf_gb  ::Array{FT,2} = zeros(FT,(325,20)) .+ 2
    "Diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    leaf_gc  ::Array{FT,2} = zeros(FT,(325,20))
    "Diffusive conductance to H₂O `[mol m⁻² s⁻¹]`"
    leaf_gh  ::Array{FT,2} = zeros(FT,(325,20))
    "Mesophyll conductance to CO₂ `[mol m⁻² s⁻¹]`"
    leaf_gm  ::Array{FT,2} = zeros(FT,(325,20)) .+ 1
    "Maximal stomatal conductance to H₂O `[mol m⁻² s⁻¹]`"
    leaf_gmax::Array{FT,2} = zeros(FT,(325,20)) .+ 1
    "Minimal stomatal conductance to H₂O `[mol m⁻² s⁻¹]`"
    leaf_gmin::Array{FT,2} = zeros(FT,(325,20)) .+ FT(0.01)
    "Stomatal conductance to H₂O `[mol m⁻² s⁻¹]`"
    leaf_gsw ::Array{FT,2} = zeros(FT,(325,20)) .+ FT(0.01)
    "Transpiration rate `[mol s⁻¹]`"
    leaf_q   ::Array{FT,2} = zeros(FT,(325,20))

    # Leaf gas exchange --- photosynthesis
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    leaf_ac   ::Array{FT,2} = zeros(FT,(325,20))
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    leaf_ag   ::Array{FT,2} = zeros(FT,(325,20))
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    leaf_aj   ::Array{FT,2} = zeros(FT,(325,20))
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    leaf_an   ::Array{FT,2} = zeros(FT,(325,20))
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    leaf_ap   ::Array{FT,2} = zeros(FT,(325,20))
    "Electron transport `[μmol m⁻² s⁻¹]`"
    leaf_j    ::Array{FT,2} = zeros(FT,(325,20))
    "Maximal electron transport rate `[μmol m⁻² s⁻¹]`"
    leaf_jm   ::Array{FT,2} = zeros(FT,(325,20))
    "Maximal electron transport rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    leaf_jm25 ::Array{FT,2} = zeros(FT,(325,20))
    "RubisCO coefficient Kc `[Pa]`"
    leaf_Kc   ::Array{FT,2} = zeros(FT,(325,20))
    "Michaelis-Menten's coefficient `[Pa]`"
    leaf_Km   ::Array{FT,2} = zeros(FT,(325,20))
    "RubisCO coefficient Ko `[Pa]`"
    leaf_Ko   ::Array{FT,2} = zeros(FT,(325,20))
    "PEP coefficient Ko `[Pa]`"
    leaf_kpep ::Array{FT,2} = zeros(FT,(325,20))
    "APAR for each leaf `[μmol m⁻² s⁻¹]`"
    leaf_apar ::Array{FT,2} = zeros(FT,(325,20)) .+ 100
    "Internal CO₂ partial pressure `[Pa]`"
    leaf_pi   ::Array{FT,2} = zeros(FT,(325,20))
    "Respiration rate `[μmol m⁻² s⁻¹]`"
    leaf_rd   ::Array{FT,2} = zeros(FT,(325,20))
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    leaf_rd25 ::Array{FT,2} = zeros(FT,(325,20))
    "Maximal carboxylation rate `[μmol m⁻² s⁻¹]`"
    leaf_vcm  ::Array{FT,2} = zeros(FT,(325,20))
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    leaf_vcm25::Array{FT,2} = zeros(FT,(325,20))
    "Maximal PEP carboxylation rate `[μmol m⁻² s⁻¹]`"
    leaf_vpm  ::Array{FT,2} = zeros(FT,(325,20))
    "Maximal PEP carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    leaf_vpm25::Array{FT,2} = zeros(FT,(325,20))
    "CO₂ compensation point with the absence of Rd `[Pa]`"
    leaf_Γs   ::Array{FT,2} = zeros(FT,(325,20))

    # Leaf fluorescence
end




function plant_initialization_canopy!(plant::PlantParaSet{FT}, lai::FT) where {FT<:AbstractFloat}
    @unpack ga,n_canopy,n_leaf = plant

    _la_layer = lai * ga / n_canopy
    for i in n_canopy
        _la_sun   = _la_layer * i / (n_canopy+1)
        _la_shade = _la_layer - _la_sun
        plant.leaf_la = FT.([ones(n_leaf-1)/(n_leaf-1)*_la_sun; _la_shade])
    end
end
=#