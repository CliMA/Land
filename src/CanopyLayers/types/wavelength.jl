###############################################################################
#
# Wave length parameter set
#
###############################################################################
"""
    mutable struct WaveLengths{FT}

Struct for pre-set wave length parameters

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct WaveLengths{FT}
    opti_file::String = LAND_2021

    # Wave length (WL) boundaries
    "Minimal WL for PAR `[nm]`"
    minwlPAR::FT = 300
    "Maximal WL for PAR `[nm]`"
    maxwlPAR::FT = 750
    "Minimal WL for NIR `[nm]`"
    minwlNIR::FT = 700
    "Maximal WL for NIR `[nm]`"
    maxwlNIR::FT = 2500
    "Minimal WL for SIF excitation `[nm]`"
    minwle::FT = 300
    "Maximal WL for SIF excitation `[nm]`"
    maxwle::FT = 750
    "Minimal WL for SIF emission/fluorescence `[nm]`"
    minwlf::FT = 640
    "Maximal WL for SIF emission/fluorescence `[nm]` "
    maxwlf::FT = 850

    # Wave length lists
    "Differential wavelength"
    dWL::Vector{FT} = read_nc(opti_file, "WL_UPPER") .- read_nc(opti_file, "WL_LOWER")

    "Leaf optical parameter set"
    optis::LeafOpticals = LeafOpticals{FT}(opti_file = opti_file)
    "Wave length `[nm]`"
    WL::Vector{FT} = read_nc(opti_file, "WL")

    "Index of WLE in WL"
    iWLE::Vector{Int} = findall( (WL .>= minwle) .& (WL .<= maxwle) )
    "Index of WLF in WL"
    iWLF::Vector{Int} = findall( (WL .>= minwlf) .& (WL .<= maxwlf) )
    "index of wlPAR in WL"
    iPAR::Vector{Int} = findall( (WL .>= minwlPAR) .& (WL .<= maxwlPAR) )
    "index of wlPAR in WL for 700 nm (regular definition)"
    iPAR_700::Vector{Int} = findall( (WL .>= minwlPAR) .& (WL .<= 700) )
    "index of wlNIR in WL"
    iNIR::Vector{Int} = findall( (WL .>= minwlNIR) .& (WL .<= maxwlNIR) )

    "excitation wave length `[nm]`"
    WLE::Vector{FT} = WL[iWLE]
    "Fluorescence wave length `[nm]`"
    WLF::Vector{FT} = WL[iWLF]
    "Wave length for PAR"
    WL_iPAR::Vector{FT} = WL[iPAR];
    "Differential wave length for PAR"
    dWL_iPAR::Vector{FT} = dWL[iPAR];
    "Differential wave length for PAR"
    dWL_iPAR_700::Vector{FT} = dWL[iPAR_700];
    "Differential wave length for iWLE"
    dWL_iWLE::Vector{FT} = dWL[iWLE];

    # local storage of dimension information
    "Length of WL_iPAR"
    nPAR::Int = length(iPAR)
    "Length of WL"
    nWL::Int = length(WL)
    "length of WLE"
    nWLE::Int = length(iWLE)
    "length of WLF"
    nWLF::Int = length(iWLF)
end
