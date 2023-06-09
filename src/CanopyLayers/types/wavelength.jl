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
    # Wave length (WL) boundaries
    "Minimal WL for PAR `[nm]`"
    minwlPAR::FT = FT(400.0)
    "Maximal WL for PAR `[nm]`"
    maxwlPAR::FT = FT(700.0)
    "Minimal WL for NIR `[nm]`"
    minwlNIR::FT = FT(700.0)
    "Maximal WL for NIR `[nm]`"
    maxwlNIR::FT = FT(2500.0)
    "Minimal WL for SIF excitation `[nm]`"
    minwle::FT = FT(400.0)
    "Maximal WL for SIF excitation `[nm]`"
    maxwle::FT = FT(750.0)
    "Minimal WL for SIF emission/fluorescence `[nm]`"
    minwlf::FT = FT(640.0)
    "Maximal WL for SIF emission/fluorescence `[nm]` "
    maxwlf::FT = FT(850.0)

    # Wave length lists
    "Standard wave length `[nm]`"
    sWL::Vector{FT} = [collect(FT(400.0):FT(10.0):FT( 650.1)); collect(FT(655.0):FT( 5.0):FT( 770.1)); collect(FT(780.0):FT(25.0):FT(2400.1))]
    "Differential wavelength"
    dWL::Vector{FT} = diff(sWL)

    "Leaf optical parameter set"
    optis::LeafOpticals = LeafOpticals{FT}()
    "Wave length `[nm]`"
    WL::Vector{FT}  = optis.lambda

    "Index of WLE in WL"
    iWLE::Vector{Int} = findall( (WL .>= minwle) .& (WL .<= maxwle) )
    "Index of WLF in WL"
    iWLF::Vector{Int} = findall( (WL .>= minwlf) .& (WL .<= maxwlf) )
    "index of wlPAR in WL"
    iPAR::Vector{Int} = findall( (WL .>= minwlPAR) .& (WL .<= maxwlPAR) )
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
