#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-15: add function to make sure the soil layers do not over saturate or drain
#     2022-Jun-18: move function to SoilPlantAirContinuum.jl
#
#######################################################################################################################################################################################################
"""

    adjusted_time(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, δt::FT) where {FT<:AbstractFloat}

Return adjusted time that soil does not over saturate or drain, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` SPAC
- `δt` Time step

"""
function adjusted_time(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, δt::FT) where {FT<:AbstractFloat}
    @unpack SOIL = spac;

    # make sure each layer does not over-saturate or drain
    _δt = δt;
    for _i in 1:SOIL.N_LAYER
        # if top soil is saturated and there is rain, _δt will not change (the rain will be counted as runoff)
        if (SOIL.LAYERS[_i].∂θ∂t > 0) && (SOIL.LAYERS[_i].θ < SOIL.LAYERS[_i].VC.Θ_SAT)
            _δt_sat = (SOIL.LAYERS[_i].VC.Θ_SAT - SOIL.LAYERS[_i].θ) / SOIL.LAYERS[_i].∂θ∂t;
            _δt = min(_δt_sat, _δt);
        elseif SOIL.LAYERS[_i].∂θ∂t < 0
            _δt_dra = (SOIL.LAYERS[_i].VC.Θ_RES - SOIL.LAYERS[_i].θ) / SOIL.LAYERS[_i].∂θ∂t;
            _δt = min(_δt_dra, _δt);
        end;
    end;

    return _δt
end
