#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-30: add abstract type for stomatal models
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractBetaFactor:
"""
abstract type AbstractBetaFactor{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-30: add abstract type for stomatal models
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractStomataModel:
"""
abstract type AbstractStomataModel{FT<:AbstractFloat} end



mutable struct BallBerrySM{FT} <: AbstractStomataModel{FT}
    # parameters that do not change with time
    "Minimal stomatal conductance `[mol m⁻² s⁻¹]`"
    G0::FT
    "Slope of conductance-photosynthesis correlation `[-]`"
    G1::FT
    "Beta function to force stomatal response tp soil moisture"
    Β
end
