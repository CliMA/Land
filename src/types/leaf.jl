###############################################################################
#
# Leaf parameters container
#
###############################################################################
"""
    mutable struct Leaf{FT}

Struct to store leaf information.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct Leaf{FT<:AbstractFloat}
    # Temperature related
    "Temperature `[K]`"
    T::FT = FT(298.15)

    # Photosynthesis system
    "Rate constant for thermal dissipation"
    Kd       ::FT = FT(0.85)
    "Rate constant for fluorescence (const)"
    Kf       ::FT = FT(0.05)
    "NPQ rate constant (initially zero)"
    Kn       ::FT = FT(0)
    "Rate constant for photochemistry (all reaction centers open)"
    Kp       ::FT = FT(4.0)
    "max PSII yield (Kn=0, all RC open)"
    maxPSII  ::FT = Kp/(Kp+Kf+Kd)
    "Fraction of absorbed light used by PSII ETR"
    PSII_frac::FT = FT(0.5)

    # CO₂ pressures
    "Leaf internal CO₂ partial pressure `[Pa]`"
    p_i  ::FT = FT(10)
    "Leaf surface CO₂ partial pressure `[Pa]`"
    p_s  ::FT = FT(40)
    "Saturation H₂O vapor pressure `[Pa]`"
    p_sat::FT = saturation_vapor_pressure(T)
    "Leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_bc ::FT = FT(3/1.35)
    "Leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_lc ::FT = FT(0.01)

    # Photosynthesis related
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ac     ::FT = FT(0)
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Aj     ::FT = FT(0)
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ag     ::FT = FT(0)
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    An     ::FT = FT(0)
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ap     ::FT = FT(0)
    "Electron transport `[μmol m⁻² s⁻¹]`"
    J      ::FT = FT(0)
    "Potential Electron Transport Rate `[μmol m⁻² s⁻¹]`"
    J_pot  ::FT = FT(0)
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

    # Fluorescence related
    "Total efficiency, incl. photorespiration `[mol CO₂ mol⁻¹ e-]`"
    CO₂_per_electron::FT = FT(1/6)
    "dark adapted yield (`Kp=0`)"
    Fm              ::FT = FT(0)
    "light adapted yield (`Kp=0`)"
    Fm′             ::FT = FT(0)
    "dark-adapted fluorescence yield (`Kp=max`)"
    Fo              ::FT = FT(0)
    "light-adapted fluorescence yield in the dark (`Kp=max`)"
    Fo′             ::FT = FT(0)
    "Actual electron transport rate `[μmol m⁻² s⁻¹]`"
    Ja              ::FT = FT(0)
    "Non-Photochemical quenching "
    NPQ             ::FT = FT(0)
    "Photochemical quenching"
    qQ              ::FT = FT(0)
    "energy quenching"
    qE              ::FT = FT(0)
    "PSII yield"
    φ               ::FT = FT(0)
    "Steady-state (light-adapted) yield (aka Fs)"
    ϕs              ::FT = FT(0)

    # Environment related
    "Absorbed photosynthetic active radiation `[μmol m⁻² s⁻¹]`"
    APAR::FT = FT(100)
end
