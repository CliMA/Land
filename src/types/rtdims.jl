###############################################################################
#
# RT Dimensions
#
###############################################################################
"""
    mutable struct RTDimentions

Struct for leaf optical properties using Array

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct RTDimentions
    "Number of azimuth angles"
    nAzi  ::Int = 36
    "Number of inclination agles"
    nIncl ::Int = 9
    "Number of canopy layers"
    nLayer::Int = 5
    "Number of canopy layer boundaries nLayer+1"
    nLevel::Int = nLayer+1
    "Number of PAR wave lengths"
    nPAR  ::Int = 10
    "Number of wave lengths"
    nWL   ::Int = 10
    "Number of wave length for excitation"
    nWLE  ::Int = 10
    "Number of wave lengths for SIF"
    nWLF  ::Int = 10
end
