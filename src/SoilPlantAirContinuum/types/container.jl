###############################################################################
#
# Containers for Simulation results
#
###############################################################################
"""
    mutable struct SPACContainer1L{FT}

Struct that contains 1-layer gas exchange information.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct SPACContainer1L{FT<:AbstractFloat}
    "Mean gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    ag::FT = FT(0)
    "Mean net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    an::FT = FT(0)
    "Leaf internal CO₂ partial pressure `[Pa]`"
    c ::FT = FT(0)
    "Flow rate per basal area `[mol s⁻¹ m⁻²]`"
    e ::FT = FT(0)
    "Leaf diffusive conductance to H₂O `[mol m⁻² s⁻¹]`"
    gh::FT = FT(0)
    "Xylem end pressure `[MPa]`"
    p ::FT = FT(0)
    "Leaf temperature `[K]`"
    t ::FT = FT(0)
end




"""
    mutable struct SPACContainer2L{FT}

Struct that contains 2-layer gas exchange information.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct SPACContainer2L{FT<:AbstractFloat}
    # gas exchange information
    cont_sh::SPACContainer1L{FT} = SPACContainer1L{FT}();
    cont_sl::SPACContainer1L{FT} = SPACContainer1L{FT}();

    # partition information
    "Shaded layer fraction"
    frac_sh::FT = FT(0)
    "Sunlit layer fraction"
    frac_sl::FT = FT(0)
    "Shaded layer leaf area `[m²]`"
    la_sh  ::FT = FT(0)
    "Sunlit layer leaf area `[m²]`"
    la_sl  ::FT = FT(0)
    "Shaded layer LAI"
    lai_sh ::FT = FT(0)
    "Sunlit layer LAI"
    lai_sl ::FT = FT(0)
    "Shaded layer PAR `[μmol m⁻² s⁻¹]`"
    par_sh ::FT = FT(0)
    "Sunlit layer PAR `[μmol m⁻² s⁻¹]`"
    par_sl ::FT = FT(0)
    "Shaded layer absorbed energy `[W m⁻²]`"
    rad_sh ::FT = FT(0)
    "Sunlit layer absorbed energy `[W m⁻²]`"
    rad_sl ::FT = FT(0)
end
