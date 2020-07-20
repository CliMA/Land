###############################################################################
#
# Wave length parameter set
#
###############################################################################
"""
    struct WaveLengths{FT}

Struct for pre-set wave length parameters.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct WaveLengths{FT}
    # Wave length (WL) boundaries
    "Minimal WL for PAR `[nm]`"
    minwlPAR::FT = FT(400.0)
    "Maximal WL for PAR `[nm]`"
    maxwlPAR::FT = FT(700.0)
    "Minimal WL for SIF excitation `[nm]`"
    minwle  ::FT = FT(400.0)
    "Maximal WL for SIF excitation `[nm]`"
    maxwle  ::FT = FT(750.0)
    "Minimal WL for SIF emission/fluorescence `[nm]`"
    minwlf  ::FT = FT(640.0)
    "Maximal WL for SIF emission/fluorescence `[nm]` "
    maxwlf  ::FT = FT(850.0)

    # Wave length lists
    "Standard wave length `[nm]`"
    swl::Array{FT,1} = [collect(FT(400.0):FT(10.0):FT( 650.1));
                        collect(FT(655.0):FT( 5.0):FT( 770.1));
                        collect(FT(780.0):FT(25.0):FT(2400.1))]
    "Differential wavelength"
    dwl::Array{FT,1} = diff(swl)

    "Leaf optical parameter set"
    optis::LeafOpticals = create_leaf_opticals(swl, file_Opti)

    "Wave length `[nm]`"
    wl::Array{FT,1} = optis.lambda

    "Length of wl"
    nwl ::Int         = length(wl)
    "Index of wle in wl"
    Iwle::Array       = findall( (wl .>= minwle) .& (wl .<= maxwle) )
    "Index of wlf in wl"
    Iwlf::Array       = findall( (wl .>= minwlf) .& (wl .<= maxwlf) )
    "length of wle"
    nWlE::Int         = length(Iwle)
    "length of wlf"
    nWlF::Int         = length(Iwlf)
    "index of wlPAR in wl"
    iPAR::Array       = findall( (wl .>= minwlPAR) .& (wl .<= maxwlPAR) )
    "excitation wave length `[nm]`"
    wle ::Array{FT,1} = wl[Iwle]
    "Fluorescence wave length `[nm]`"
    wlf ::Array{FT,1} = wl[Iwlf]
end
