###############################################################################
#
# CanopyLayer parameters container
#
###############################################################################
"""
    struct CanopyLayer{FT}

Struct to store leaf information (multi-dimensional).

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct CanopyLayer{FT<:AbstractFloat}
    # CanopyLayer photosynthesis system
    "leaf photosynthesis system"
    ps  ::Leaf{FT} = Leaf{FT}()
    "Memory leaf photosynthesis system"
    ps_m::Leaf{FT} = Leaf{FT}()
    "Total leaf area `[m²]`"
    LA    ::FT = 150
    "Leaf area index in the layer"
    LAI   ::FT = 3
    "Total leaf area index in the layer"
    tLAI  ::FT = 3
    "Memory APAR `[μmol m⁻² s⁻¹]`"
    APAR_m::FT = 500
    "Memory environment"
    envir_m::AirLayer{FT} = AirLayer{FT}()

    # Number of leaves per canopy layer
    n_leaf::Int = 325

    # Temperature related, different for each leaf
    "Sensible Heat Flux `[W m⁻²]`"
    H ::Vector{FT} = zeros(FT, n_leaf)
    "Latent Heat Flux `[W m⁻²]`"
    LE::Vector{FT} = zeros(FT, n_leaf)
    "Net Radiation Balance `[W m⁻²]`"
    Rn::Vector{FT} = zeros(FT, n_leaf)

    # Tempearture related, same for all leaves
    "Latent Heat of evaporation `[J mol⁻¹]`"
    LV   ::FT = latent_heat_vapor(T₂₅(FT)) * M_H₂O(FT)
    "Temperature `[K]`"
    T    ::FT = T₂₅(FT)
    "Old temperature `[K]`"
    T_old::FT = FT(0)
    "Leaf width `[m]`"
    width::FT = FT(0.05)

    # Diffusive conductances, same for all leaves
    "Boundary layer conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_bc::Vector{FT} = zeros(FT, n_leaf) .+ FT(3/1.35)
    "Boundary layer conductance to heat `[mol m⁻² s⁻¹]`"
    g_bh::Vector{FT} = zeros(FT, n_leaf) .+ FT(3)
    "Boundary layer conductance to H₂O `[mol m⁻² s⁻¹]`"
    g_bw::Vector{FT} = zeros(FT, n_leaf) .+ FT(3)
    "Leaf diffusive conductance to water CO₂ `[mol m⁻² s⁻¹]`"
    g_lc::Vector{FT} = zeros(FT, n_leaf) .+ FT(0.08032)
    "Leaf diffusive conductance to water H₂O `[mol m⁻² s⁻¹]`"
    g_lw::Vector{FT} = zeros(FT, n_leaf) .+ FT(0.1519)
    "Mesophyll conductance for CO₂ `[mol m⁻² s⁻¹]`"
    g_m ::Vector{FT} = zeros(FT, n_leaf) .+ FT(0.5)
    "Stomatal conductance to water CO₂ `[mol m⁻² s⁻¹]`"
    g_sc::Vector{FT} = zeros(FT, n_leaf) .+ FT(0.1)
    "Stomatal conductance to water H₂O `[mol m⁻² s⁻¹]`"
    g_sw::Vector{FT} = zeros(FT, n_leaf) .+ FT(0.16)

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
    g_min  ::FT = FT(0.001)
    "Minimal leaf diffusive conductance at 298.15 K `[mol m⁻² s⁻¹]`"
    g_min25::FT = FT(0.001)

    # CO₂ and H₂O pressures, different for each leaf
    "Leaf internal CO₂ partial pressure `[Pa]`"
    p_i::Vector{FT} = zeros(FT, n_leaf) .+ FT(10)
    "Leaf surface CO₂ partial pressure `[Pa]`"
    p_s::Vector{FT} = zeros(FT, n_leaf) .+ FT(40)

    # CO₂ and H₂O pressures, same for all leaves
    "Leaf saturation vapor pressure `[Pa]`"
    p_sat::FT = saturation_vapor_pressure(T)

    # Photosynthesis related, different for each leaf
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ac   ::Vector{FT} = zeros(FT, n_leaf)
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Aj   ::Vector{FT} = zeros(FT, n_leaf)
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ag   ::Vector{FT} = zeros(FT, n_leaf)
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    An   ::Vector{FT} = zeros(FT, n_leaf)
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ap   ::Vector{FT} = zeros(FT, n_leaf)
    "Electron transport `[μmol m⁻² s⁻¹]`"
    J    ::Vector{FT} = zeros(FT, n_leaf)
    "Potential Electron Transport Rate `[μmol m⁻² s⁻¹]`"
    J_pot::Vector{FT} = zeros(FT, n_leaf)

    # Fluorescence related, different for each leaf
    "Total efficiency, incl. photorespiration `[mol CO₂ mol⁻¹ e-]`"
    e2c::Vector{FT} = zeros(FT, n_leaf) .+ FT(1/6)
    "light adapted yield (`Kp=0`)"
    Fm′::Vector{FT} = zeros(FT, n_leaf)
    "light-adapted fluorescence yield in the dark (`Kp=max`)"
    Fo′::Vector{FT} = zeros(FT, n_leaf)
    "Actual electron transport rate `[μmol m⁻² s⁻¹]`"
    Ja ::Vector{FT} = zeros(FT, n_leaf)
    "Non-Photochemical quenching "
    NPQ::Vector{FT} = zeros(FT, n_leaf)
    "Photochemical quenching"
    qQ ::Vector{FT} = zeros(FT, n_leaf)
    "energy quenching"
    qE ::Vector{FT} = zeros(FT, n_leaf)
    "PSII yield"
    φ  ::Vector{FT} = zeros(FT, n_leaf)
    "Steady-state (light-adapted) yield (aka Fs)"
    φs ::Vector{FT} = zeros(FT, n_leaf)

    # Fluorescence related, same for all leaves
    "dark adapted yield (`Kp=0`)"
    Fm::FT = FT(0)
    "dark-adapted fluorescence yield (`Kp=max`)"
    Fo::FT = FT(0)

    # Environment related, different for each leaf
    "Absorbed photosynthetic active radiation `[μmol m⁻² s⁻¹]`"
    APAR::Vector{FT} = zeros(FT, n_leaf) .+ 100
    "Leaf area fractions"
    LAIx::Vector{FT} = ones(FT, n_leaf) ./ n_leaf;

    # Stomtal optimization related, different for each leaf
    "Maximal photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_max::Vector{FT} = zeros(FT, n_leaf)
    "Flow rate `[mol m⁻² s⁻¹]`"
    e    ::Vector{FT} = zeros(FT, n_leaf)

    # Stomtal optimization related, same for all leaves
    "Critical flow rate `[mol m⁻² s⁻¹]`"
    ec    ::FT = FT(2e-9)
    "Maximal hydraulic conductance ratio"
    kr_max::FT = FT(1)
    "Base xylem pressre `[MPa]`"
    p_ups ::FT = FT(0)
    "Base xylem pressre memory `[MPa]`"
    p_old ::FT = FT(1)
    "Fitness factor for nighttime stomtal conductance"
    ff    ::FT = FT(0.15)

    # time constant related to prognostic stomatal conductance
    "τ for empirical stomatal models `[-]`, Δg/Δt = (g_ss - gsw) / τ"
    τ_esm::FT = 600
    "τ for optimal stomatal models `[μmol⁻¹]`, Δg/Δt = (∂A/∂E - ∂Θ/∂E) * τ"
    τ_osm::FT = FT(1e-6)
    "τ for nighttime optimal stomatal models `[μmol⁻¹]`"
    τ_noc::FT = FT(2e-6)
end
