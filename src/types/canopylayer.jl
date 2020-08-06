###############################################################################
#
# CanopyLayer parameters container
#
###############################################################################
"""
    struct CanopyLayer{FT}

Struct to store leaf information (multi-dimensional).

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct CanopyLayer{FT<:AbstractFloat}
    # CanopyLayer photosynthesis system
    "leaf photosynthesis system"
    ps ::Leaf{FT} = Leaf{FT}()
    "Total leaf area `[m²]`"
    LA ::FT = FT(150)
    "Total leaf area index"
    LAI::FT = FT(3)

    # Number of leaves per canopy layer
    n_leaf::Int = 325

    # Temperature related, different for each leaf
    "Sensible Heat Flux `[W m⁻²]`"
    H ::Array{FT,1} = zeros(FT, n_leaf)
    "Latent Heat Flux `[W m⁻²]`"
    LE::Array{FT,1} = zeros(FT, n_leaf)
    "Net Radiation Balance `[W m⁻²]`"
    Rn::Array{FT,1} = zeros(FT, n_leaf)

    # Tempearture related, same for all leaves
    "Latent Heat of evaporation `[J mol⁻¹]`"
    LV   ::FT = latent_heat_vapor(K_25(FT)) * 1000 / MOLMASS_WATER(FT)
    "Temperature `[K]`"
    T    ::FT = K_25(FT)
    "Old temperature `[K]`"
    T_old::FT = FT(0)
    "Leaf width `[m]`"
    width::FT = FT(0.05)

    # Photosynthesis system, different for each leaf
    "NPQ rate constant (initially zero)"
    Kn     ::Array{FT,1} = zeros(FT, n_leaf)
    "Rate constant for photochemistry (all reaction centers open)"
    Kp     ::Array{FT,1} = zeros(FT, n_leaf) .+ 4

    # Photosynthesis system, same for all leaves
    "Rate constant for thermal dissipation"
    Kd       ::FT = FT(0.85)
    "Rate constant for fluorescence (const)"
    Kf       ::FT = FT(0.05)
    "Maximal rate constant for photochemistry (all reaction centers open)"
    Kpmax    ::FT = FT(4)
    "max PSII yield (Kn=0, all RC open)"
    maxPSII  ::FT = Kpmax / (Kpmax + Kf +Kd)
    "Fraction of absorbed light used by PSII ETR"
    PSII_frac::FT = FT(0.5)

    # Diffusive conductances, same for all leaves
    "Boundary layer conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_bc::Array{FT,1} = zeros(FT, n_leaf) .+ FT(3/1.35)
    "Boundary layer conductance to heat `[mol m⁻² s⁻¹]`"
    g_bh::Array{FT,1} = zeros(FT, n_leaf) .+ FT(3)
    "Boundary layer conductance to H₂O `[mol m⁻² s⁻¹]`"
    g_bw::Array{FT,1} = zeros(FT, n_leaf) .+ FT(3)
    "Leaf diffusive conductance to water CO₂ `[mol m⁻² s⁻¹]`"
    g_lc::Array{FT,1} = zeros(FT, n_leaf) .+ FT(0.08032)
    "Leaf diffusive conductance to water H₂O `[mol m⁻² s⁻¹]`"
    g_lw::Array{FT,1} = zeros(FT, n_leaf) .+ FT(0.1519)
    "Mesophyll conductance for CO₂ `[mol m⁻² s⁻¹]`"
    g_m ::Array{FT,1} = zeros(FT, n_leaf) .+ FT(0.5)
    "Stomatal conductance to water CO₂ `[mol m⁻² s⁻¹]`"
    g_sc::Array{FT,1} = zeros(FT, n_leaf) .+ FT(0.1)
    "Stomatal conductance to water H₂O `[mol m⁻² s⁻¹]`"
    g_sw::Array{FT,1} = zeros(FT, n_leaf) .+ FT(0.16)

    # Diffusive conductances, same for all leaves
    "Gias correction constant"
    g_ias_c::FT = FT(0)
    "Gias correction exponent"
    g_ias_e::FT = FT(0.3)
    "Maximal leaf diffusive conductance `[mol m⁻² s⁻¹]`"
    g_max  ::FT = FT(0.8)
    "Maximal leaf diffusive conductance at 298.15 K `[mol m⁻² s⁻¹]`"
    g_max25::FT = FT(0.8)
    "Minimal leaf diffusive conductance `[mol m⁻² s⁻¹]`"
    g_min  ::FT = FT(0.025)
    "Minimal leaf diffusive conductance at 298.15 K `[mol m⁻² s⁻¹]`"
    g_min25::FT = FT(0.025)

    # CO₂ and H₂O pressures, different for each leaf
    "Leaf internal CO₂ partial pressure `[Pa]`"
    p_i::Array{FT,1} = zeros(FT, n_leaf) .+ FT(10)
    "Leaf surface CO₂ partial pressure `[Pa]`"
    p_s::Array{FT,1} = zeros(FT, n_leaf) .+ FT(40)

    # CO₂ and H₂O pressures, same for all leaves
    "Leaf saturation vapor pressure `[Pa]`"
    p_sat::FT = saturation_vapor_pressure(T)

    # Photosynthesis related, different for each leaf
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ac   ::Array{FT,1} = zeros(FT, n_leaf)
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Aj   ::Array{FT,1} = zeros(FT, n_leaf)
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ag   ::Array{FT,1} = zeros(FT, n_leaf)
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    An   ::Array{FT,1} = zeros(FT, n_leaf)
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ap   ::Array{FT,1} = zeros(FT, n_leaf)
    "Electron transport `[μmol m⁻² s⁻¹]`"
    J    ::Array{FT,1} = zeros(FT, n_leaf)
    "Potential Electron Transport Rate `[μmol m⁻² s⁻¹]`"
    J_pot::Array{FT,1} = zeros(FT, n_leaf)

    # Photosynthesis related, same for all leaves
    "Maximal electron transport rate `[μmol m⁻² s⁻¹]`"
    Jmax   ::FT = FT(120)
    "Maximal electron transport rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Jmax25 ::FT = FT(120)
    "RubisCO coefficient Kc `[Pa]`"
    Kc     ::FT = FT(0)
    "RubisCO coefficient Ko `[Pa]`"
    Ko     ::FT = FT(0)
    "PEP coefficient Ko `[Pa]`"
    Kpep   ::FT = FT(0)
    "Michaelis-Menten's coefficient `[Pa]`"
    Km     ::FT = FT(0)
    "Respiration rate `[μmol m⁻² s⁻¹]`"
    Rd     ::FT = FT(1)
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Rd25   ::FT = FT(1)
    "Maximal carboxylation rate `[μmol m⁻² s⁻¹]`"
    Vcmax  ::FT = FT(60)
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Vcmax25::FT = FT(60)
    "Maximal PEP carboxylation rate `[μmol m⁻² s⁻¹]`"
    Vpmax  ::FT = FT(120)
    "Maximal PEP carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Vpmax25::FT = FT(120)
    "CO₂ compensation point with the absence of Rd `[Pa]`"
    Γ_star ::FT = FT(0)

    # Fluorescence related, different for each leaf
    "Total efficiency, incl. photorespiration `[mol CO₂ mol⁻¹ e-]`"
    CO₂_per_electron::Array{FT,1} = zeros(FT, n_leaf) .+ FT(1/6)
    "light adapted yield (`Kp=0`)"
    Fm′             ::Array{FT,1} = zeros(FT, n_leaf)
    "light-adapted fluorescence yield in the dark (`Kp=max`)"
    Fo′             ::Array{FT,1} = zeros(FT, n_leaf)
    "Actual electron transport rate `[μmol m⁻² s⁻¹]`"
    Ja              ::Array{FT,1} = zeros(FT, n_leaf)
    "Non-Photochemical quenching "
    NPQ             ::Array{FT,1} = zeros(FT, n_leaf)
    "Photochemical quenching"
    qQ              ::Array{FT,1} = zeros(FT, n_leaf)
    "energy quenching"
    qE              ::Array{FT,1} = zeros(FT, n_leaf)
    "PSII yield"
    φ               ::Array{FT,1} = zeros(FT, n_leaf)
    "Steady-state (light-adapted) yield (aka Fs)"
    ϕs              ::Array{FT,1} = zeros(FT, n_leaf)

    # Fluorescence related, same for all leaves
    "dark adapted yield (`Kp=0`)"
    Fm::FT = FT(0)
    "dark-adapted fluorescence yield (`Kp=max`)"
    Fo::FT = FT(0)

    # Environment related, different for each leaf
    "Absorbed photosynthetic active radiation `[μmol m⁻² s⁻¹]`"
    APAR::Array{FT,1} = zeros(FT, n_leaf) .+ 100
    "Leaf area fractions"
    LAIx::Array{FT,1} = ones(FT, n_leaf) ./ n_leaf;

    # Stomtal optimization related, different for each leaf
    "Maximal photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_max::Array{FT,1} = zeros(FT, n_leaf)
    "Flow rate `[mol m⁻² s⁻¹]`"
    e    ::Array{FT,1} = zeros(FT, n_leaf)

    # Stomtal optimization related, same for all leaves
    "Critical flow rate `[mol m⁻² s⁻¹]`"
    ec    ::FT = FT(2e-9)
    "Maximal hydraulic conductance ratio"
    kr_max::FT = FT(1)
    "Base xylem pressre `[MPa]`"
    p_ups ::FT = FT(0)
    "Base xylem pressre memory `[MPa]`"
    p_old ::FT = FT(1)
end
