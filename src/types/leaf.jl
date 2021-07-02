###############################################################################
#
# Leaf parameters container
#
###############################################################################
"""
    mutable struct Leaf{FT}

Struct to store leaf information.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Leaf{FT<:AbstractFloat}
    # Temperature related
    "Temperature `[K]`"
    T::FT = 298.15
    "Old Temperature `[K]`, if not T, run leaf_temperature_dependence!"
    T_old::FT = 0

    # Photosynthesis system
    "Rate constant for thermal dissipation"
    Kd::FT = 0.85
    "Rate constant for fluorescence (const)"
    Kf::FT = 0.05
    "Reversible NPQ rate constant (initially zero)"
    Kr::FT = 0
    "Sustained NPQ rate constant (for seasonal changes, default is zero)"
    Ks::FT = 0
    "Rate constant for photochemistry (all reaction centers open)"
    Kp::FT = 4
    "Maximal Kp"
    Kp_max::FT = 4
    "max PSII yield (Kr=0, all RC open)"
    maxPSII::FT = Kp/(Kp+Kf+Kd)
    "Fraction of absorbed light used by PSII ETR"
    PSII_frac::FT = 0.5

    # CO₂ pressures
    "Leaf internal CO₂ partial pressure `[Pa]`"
    p_i::FT = 10
    "Leaf surface CO₂ partial pressure `[Pa]`"
    p_s::FT = 40
    "Saturation H₂O vapor pressure `[Pa]`"
    p_sat::FT = saturation_vapor_pressure(T)
    "Leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_bc::FT = 3 / 1.35
    "Leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_lc::FT = 0.01

    # Photosynthesis related
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ac::FT = 0
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Aj::FT = 0
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ag::FT = 0
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    An::FT = 0
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ap::FT = 0
    "Electron transport `[μmol m⁻² s⁻¹]`"
    J::FT = 0
    "Potential Electron Transport Rate `[μmol m⁻² s⁻¹]`"
    J_pot::FT = 0
    "Maximal electron transport rate `[μmol m⁻² s⁻¹]`"
    Jmax::FT = 120
    "Maximal electron transport rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Jmax25::FT = 120
    "RubisCO coefficient Kc `[Pa]`"
    Kc::FT = 0
    "RubisCO coefficient Ko `[Pa]`"
    Ko::FT = 0
    "PEP coefficient Ko `[Pa]`"
    Kpep::FT = 0
    "Michaelis-Menten's coefficient `[Pa]`"
    Km::FT = 0
    "Respiration rate `[μmol m⁻² s⁻¹]`"
    Rd::FT = 1
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Rd25::FT = 1
    "Maximal carboxylation rate `[μmol m⁻² s⁻¹]`"
    Vcmax::FT = 60
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Vcmax25::FT = 60
    "Maximal PEP carboxylation rate `[μmol m⁻² s⁻¹]`"
    Vpmax::FT = 120
    "Maximal PEP carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Vpmax25::FT = 120
    "CO₂ compensation point with the absence of Rd `[Pa]`"
    Γ_star::FT = 0

    # Cytochrome related
    "Total concentration of Cytochrome b₆f `[μmol m⁻²]`"
    C_b₆f::FT = 350 / 300
    "Maximal turnover rate of Cytochrome b₆f `[e⁻ s⁻¹]`"
    k_q::FT = 300
    "Maximal Cytochrome b₆f activity `[μmol e⁻ m⁻² s⁻¹]`"
    Vqmax::FT = C_b₆f * k_q
    "Rate constant of consititutive heat loss from the antennae `[s⁻¹]`"
    K_D1::FT = 5.5e8
    "rate constant of fluorescence `[s⁻¹]`"
    K_F1::FT = 5e7
    "Rate constant of regulated heat loss for PS I `[s⁻¹]`"
    K_N1::FT = 1.45e10
    "Rate constant of photochemistry for PS I `[s⁻¹]`"
    K_P1::FT = 1.45e10
    "Rate constant of photochemistry for PS II `[s⁻¹]`"
    K_P2::FT = 4.5e9
    "Rate constant of excitation sharing for PS II `[s⁻¹]`"
    K_U2::FT = 0
    "PPFD absorbed by PS I per incident PPFD"
    α_1::FT = 0.41
    "PPFD absorbed by PS II per incident PPFD"
    α_2::FT = 0.44
    "Weighting factor for PS I"
    ϵ_1::FT = 0
    "Weighting factor for PS II"
    ϵ_2::FT = 1
    "Maximal PS I photochemical yield"
    φ_P1_max::FT = K_P1 / (K_P1 + K_D1 + K_F1)
    "Coupling efficiency of cyclic electron flow `[mol ATP mol⁻¹ e⁻]`"
    n_C::FT = 1
    "Coupling efficiency of linear electron flow `[mol ATP mol⁻¹ e⁻]`"
    n_L::FT = 0.75
    "ratio between J_P700 and J_P680"
    η::FT = 0

    # related to fluorescence
    "PS II electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    J_P680_a::FT = 0
    "Rubisco limited PS II electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    J_P680_c::FT = 0
    "Light limited PS II electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    J_P680_j::FT = 0
    "Product limited PS II electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    J_P680_p::FT = 0
    "PS I electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    J_P700_a::FT = 0
    "Rubisco limited PS I electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    J_P700_c::FT = 0
    "Light limited PS I electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    J_P700_j::FT = 0
    "Product limited PS I electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    J_P700_p::FT = 0

    # Well watered condition values (to use with β function over PS)
    "Well watered maximal electron transport rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Jmax25WW::FT = Jmax25
    "Well watered respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Rd25WW::FT = Rd25
    "Well watered maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Vcmax25WW::FT = Vcmax25
    "Well watered maximal PEP carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Vpmax25WW::FT = Vpmax25

    # Fluorescence related
    "Total efficiency, incl. photorespiration `[mol CO₂ mol⁻¹ e-]`"
    e2c::FT = 1 / 6
    "dark adapted yield (`Kp=0`)"
    Fm::FT = 0
    "light adapted yield (`Kp=0`)"
    Fm′::FT = 0
    "dark-adapted fluorescence yield (`Kp=max`)"
    Fo::FT = 0
    "light-adapted fluorescence yield in the dark (`Kp=max`)"
    Fo′::FT = 0
    "Actual electron transport rate `[μmol m⁻² s⁻¹]`"
    Ja::FT = 0
    "Non-Photochemical quenching "
    NPQ::FT = 0
    "Photochemical quenching"
    qQ::FT = 0
    "energy quenching"
    qE::FT = 0
    "PSII yield"
    φ::FT = 0
    "Steady-state (light-adapted) yield (aka Fs)"
    φs::FT = 0

    # Environment related
    "Absorbed photosynthetic active radiation `[μmol m⁻² s⁻¹]`"
    APAR::FT = 100
end
