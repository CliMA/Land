###############################################################################
#
# Method to fit soil albedo
#
###############################################################################
"""
    abstract type AbstractAlbedoFitting

Hierachy of AbstractAlbedoFitting
- [`FourBandsFittingCurve`](@ref)
- [`FourBandsFittingHybrid`](@ref)
- [`FourBandsFittingPoint`](@ref)
- [`TwoBandsFittingCurve`](@ref)
- [`TwoBandsFittingHybrid`](@ref)
- [`TwoBandsFittingPoint`](@ref)
"""
abstract type AbstractAlbedoFitting end




"""
$(TYPEDEF)

Method to use all four GSV PCA bands to fit two flat lines
"""
struct FourBandsFittingCurve <: AbstractAlbedoFitting end




"""
$(TYPEDEF)

Method to use all four GSV PCA bands to fit a mean and a flat curve
"""
struct FourBandsFittingHybrid <: AbstractAlbedoFitting end




"""
$(TYPEDEF)

Method to use all four GSV PCA bands to fit 2 mean values
"""
struct FourBandsFittingPoint <: AbstractAlbedoFitting end




"""
$(TYPEDEF)

Method to use 2 GSV PCA bands to fit two flat lines
"""
struct TwoBandsFittingCurve <: AbstractAlbedoFitting end




"""
$(TYPEDEF)

Method to use 2 GSV PCA bands to fit a mean and a flat curve
"""
struct TwoBandsFittingHybrid <: AbstractAlbedoFitting end




"""
$(TYPEDEF)

Method to use 2 GSV PCA bands to fit 2 mean values
"""
struct TwoBandsFittingPoint <: AbstractAlbedoFitting end








###############################################################################
#
# Calculate soil albedo based on color class and soil moisture
#
###############################################################################
"""
    soil_albedos(color::Int, swc::FT) where {FT<:AbstractFloat}

Calculate soil broad band albedos for PAR and NIR, given
- `color` Soil color class, 1~20
- `swc` Soil volumetric water content
"""
function soil_albedos(color::Int, swc::FT) where {FT<:AbstractFloat}
    # soil class must be from 1 to 20
    @assert 1 <= color <=20;

    # calculate albedos for PAR and NIR
    _delta   = max(0, FT(0.11) - FT(0.4) * swc);
    _alb_par::FT = max(SOIL_BNDS[color,1], SOIL_BNDS[color,3] + _delta);
    _alb_nir::FT = max(SOIL_BNDS[color,2], SOIL_BNDS[color,4] + _delta);

    return _alb_par,_alb_nir
end




function soil_albedos(
            color::Int,
            rswc::FT,
            rel::Bool
) where {FT<:AbstractFloat}
    # soil class must be from 1 to 20
    @assert 1 <= color <=20;

    # calculate albedos for PAR and NIR
    _alb_par::FT = SOIL_BNDS[color,1] * (1 - rswc) + rswc * SOIL_BNDS[color,3];
    _alb_nir::FT = SOIL_BNDS[color,2] * (1 - rswc) + rswc * SOIL_BNDS[color,4];

    return _alb_par,_alb_nir
end




"""
    fit_soil_mat!(
                soil::SoilOpticals{FT},
                wls::WaveLengths{FT},
                swc::FT,
                method::AbstractAlbedoFitting;
                clm::Bool = false
    ) where {FT<:AbstractFloat}

Fit soil albedo bands parameters, given
- `soil` [`SoilOpticals`](@ref) type structure
- `wls` [`WaveLengths`](@ref) type structure
- `swc` Soil volumetric water content
- `method` [`AbstractAlbedoFitting`](@ref) type fitting method
- `clm` If true, use CLM method, else use new method
"""
function fit_soil_mat!(
            soil::SoilOpticals{FT},
            wls::WaveLengths{FT},
            swc::FT,
            method::AbstractAlbedoFitting;
            clm::Bool = false
) where {FT<:AbstractFloat}
    @unpack iWLF = wls;

    # fit the curve only if the values mismatch
    if clm
        _ref_par,_ref_nir = soil_albedos(soil.color, swc);
    else
        _ref_par,_ref_nir = soil_albedos(soil.color, swc, true);
    end;
    if soil.ρ_PAR !== _ref_par || soil.ρ_NIR !== _ref_nir
        fit_soil_mat!(soil, wls, _ref_par, _ref_nir, method);
    end

    # update the soil information
    soil.ρ_SW_SIF .= view(soil.ρ_SW, iWLF);
    soil.ε_SW .= 1 .- soil.ρ_SW;

    return nothing
end




"""
    fit_soil_mat!(
                soil::SoilOpticals{FT},
                wls::WaveLengths{FT},
                ref_PAR::FT,
                ref_NIR::FT,
                method::FourBandsFittingCurve
    ) where {FT<:AbstractFloat}

Fit soil albedo bands parameters, given
- `soil` [`SoilOpticals`](@ref) type structure
- `wls` [`WaveLengths`](@ref) type structure
- `ref_PAR` Target albedo for PAR
- `ref_NIR` Target albedo for NIR
- `method` [`FourBandsFittingCurve`](@ref) type fitting method
"""
function fit_soil_mat!(
            soil::SoilOpticals{FT},
            wls::WaveLengths{FT},
            ref_PAR::FT,
            ref_NIR::FT,
            method::FourBandsFittingCurve
) where {FT<:AbstractFloat}
    # update soil PAR and NIR albedo
    soil.ρ_PAR = ref_PAR;
    soil.ρ_NIR = ref_NIR;

    # solve for weights using pinv
    soil.ρ_SW_raw[32:end] .= ref_NIR;
    soil.ρ_SW_raw[1:31] .= ref_PAR;
    soil.SW_vec_4 .= pinv(soil.SW_mat_raw_4) * soil.ρ_SW_raw;

    # update vectors in soil
    mul!(soil.ρ_SW, soil.SW_mat_4, soil.SW_vec_4);

    return nothing
end




"""
    fit_soil_mat!(
                soil::SoilOpticals{FT},
                wls::WaveLengths{FT},
                ref_PAR::FT,
                ref_NIR::FT,
                method::FourBandsFittingHybrid
    ) where {FT<:AbstractFloat}

Fit soil albedo bands parameters, given
- `soil` [`SoilOpticals`](@ref) type structure
- `wls` [`WaveLengths`](@ref) type structure
- `ref_PAR` Target albedo for PAR
- `ref_NIR` Target albedo for NIR
- `method` [`FourBandsFittingHybrid`](@ref) type fitting method
"""
function fit_soil_mat!(
            soil::SoilOpticals{FT},
            wls::WaveLengths{FT},
            ref_PAR::FT,
            ref_NIR::FT,
            method::FourBandsFittingHybrid
) where {FT<:AbstractFloat}
    # update soil PAR and NIR albedo
    soil.ρ_PAR = ref_PAR;
    soil.ρ_NIR = ref_NIR;

    # solve for weights using pinv
    soil.ρ_SW_raw[32:end] .= ref_NIR;
    soil.ρ_SW_raw[1:31] .= ref_PAR;
    soil.SW_vec_4 .= pinv(soil.SW_mat_raw_4) * soil.ρ_SW_raw;

    # solve for weights
    @inline _fit(x::Vector{FT}) where {FT<:AbstractFloat} = (
        mul!(soil.ρ_SW_raw, soil.SW_mat_raw_4, x);
        _diff = ( mean(soil.ρ_SW_raw[1:31]) - ref_PAR ) ^ 2 +
                mean( abs.(soil.ρ_SW_raw[32:end] .- ref_NIR) ) ^ 2;
        return -_diff
    );

    _ms = ReduceStepMethodND{FT}(x_mins = FT[-2,-2,-2,-2],
                                 x_maxs = FT[2,2,2,2],
                                 x_inis = soil.SW_vec_4,
                                 Δ_inis = FT[0.1,0.1,0.1,0.1]);
    _tol = SolutionToleranceND{FT}(FT[0.001,0.001,0.001,0.001], 50);
    _sol = find_peak(_fit, _ms, _tol);
    soil.SW_vec_4 .= _sol;

    # update vectors in soil
    mul!(soil.ρ_SW, soil.SW_mat_4, soil.SW_vec_4);

    return nothing
end




"""
    fit_soil_mat!(
                soil::SoilOpticals{FT},
                wls::WaveLengths{FT},
                ref_PAR::FT,
                ref_NIR::FT,
                method::FourBandsFittingPoint
    ) where {FT<:AbstractFloat}

Fit soil albedo bands parameters, given
- `soil` [`SoilOpticals`](@ref) type structure
- `wls` [`WaveLengths`](@ref) type structure
- `ref_PAR` Target albedo for PAR
- `ref_NIR` Target albedo for NIR
- `method` [`FourBandsFittingPoint`](@ref) type fitting method
"""
function fit_soil_mat!(
            soil::SoilOpticals{FT},
            wls::WaveLengths{FT},
            ref_PAR::FT,
            ref_NIR::FT,
            method::FourBandsFittingPoint
) where {FT<:AbstractFloat}
    # update soil PAR and NIR albedo
    soil.ρ_PAR = ref_PAR;
    soil.ρ_NIR = ref_NIR;

    # solve for weights using pinv
    soil.ρ_SW_raw[32:end] .= ref_NIR;
    soil.ρ_SW_raw[1:31] .= ref_PAR;
    soil.SW_vec_4 .= pinv(soil.SW_mat_raw_4) * soil.ρ_SW_raw;

    # solve for weights
    @inline _fit(x::Vector{FT}) where {FT<:AbstractFloat} = (
        mul!(soil.ρ_SW_raw, soil.SW_mat_raw_4, x);
        _diff = ( mean(soil.ρ_SW_raw[1:31]) - ref_PAR ) ^ 2 +
                ( mean(soil.ρ_SW_raw[32:end]) - ref_NIR ) ^ 2;
        return -_diff
    );

    _ms = ReduceStepMethodND{FT}(x_mins = FT[-2,-2,-2,-2],
                                 x_maxs = FT[2,2,2,2],
                                 x_inis = soil.SW_vec_4,
                                 Δ_inis = FT[0.1,0.1,0.1,0.1]);
    _tol = SolutionToleranceND{FT}(FT[0.001,0.001,0.001,0.001], 50);
    _sol = find_peak(_fit, _ms, _tol);
    soil.SW_vec_4 .= _sol;

    # update vectors in soil
    mul!(soil.ρ_SW, soil.SW_mat_4, soil.SW_vec_4);

    return nothing
end




"""
    fit_soil_mat!(
                soil::SoilOpticals{FT},
                wls::WaveLengths{FT},
                ref_PAR::FT,
                ref_NIR::FT,
                method::TwoBandsFittingCurve
    ) where {FT<:AbstractFloat}

Fit soil albedo bands parameters, given
- `soil` [`SoilOpticals`](@ref) type structure
- `wls` [`WaveLengths`](@ref) type structure
- `ref_PAR` Target albedo for PAR
- `ref_NIR` Target albedo for NIR
- `method` [`TwoBandsFittingCurve`](@ref) type fitting method
"""
function fit_soil_mat!(
            soil::SoilOpticals{FT},
            wls::WaveLengths{FT},
            ref_PAR::FT,
            ref_NIR::FT,
            method::TwoBandsFittingCurve
) where {FT<:AbstractFloat}
    # update soil PAR and NIR albedo
    soil.ρ_PAR = ref_PAR;
    soil.ρ_NIR = ref_NIR;

    # solve for weights using pinv
    soil.ρ_SW_raw[32:end] .= ref_NIR;
    soil.ρ_SW_raw[1:31] .= ref_PAR;
    soil.SW_vec_2 .= pinv(soil.SW_mat_raw_2) * soil.ρ_SW_raw;

    # update vectors in soil
    mul!(soil.ρ_SW, soil.SW_mat_2, soil.SW_vec_2);

    return nothing
end




"""
    fit_soil_mat!(
                soil::SoilOpticals{FT},
                wls::WaveLengths{FT},
                ref_PAR::FT,
                ref_NIR::FT,
                method::TwoBandsFittingHybrid
    ) where {FT<:AbstractFloat}

Fit soil albedo bands parameters, given
- `soil` [`SoilOpticals`](@ref) type structure
- `wls` [`WaveLengths`](@ref) type structure
- `ref_PAR` Target albedo for PAR
- `ref_NIR` Target albedo for NIR
- `method` [`TwoBandsFittingHybrid`](@ref) type fitting method
"""
function fit_soil_mat!(
            soil::SoilOpticals{FT},
            wls::WaveLengths{FT},
            ref_PAR::FT,
            ref_NIR::FT,
            method::TwoBandsFittingHybrid
) where {FT<:AbstractFloat}
    # update soil PAR and NIR albedo
    soil.ρ_PAR = ref_PAR;
    soil.ρ_NIR = ref_NIR;

    # solve for weights using pinv
    soil.ρ_SW_raw[32:end] .= ref_NIR;
    soil.ρ_SW_raw[1:31] .= ref_PAR;
    soil.SW_vec_2 .= pinv(soil.SW_mat_raw_2) * soil.ρ_SW_raw;

    # solve for weights
    @inline _fit(x::Vector{FT}) where {FT<:AbstractFloat} = (
        mul!(soil.ρ_SW_raw, soil.SW_mat_raw_2, x);
        _diff = ( mean(soil.ρ_SW_raw[1:31]) - ref_PAR ) ^ 2 +
                mean( abs.(soil.ρ_SW_raw[32:end] .- ref_NIR) ) ^ 2;
        return -_diff
    );

    _ms = ReduceStepMethodND{FT}(x_mins = FT[-2,-2],
                                 x_maxs = FT[2,2],
                                 x_inis = soil.SW_vec_2,
                                 Δ_inis = FT[0.1,0.1]);
    _tol = SolutionToleranceND{FT}(FT[0.001,0.001], 50);
    _sol = find_peak(_fit, _ms, _tol);
    soil.SW_vec_2 .= _sol;

    # update vectors in soil
    mul!(soil.ρ_SW, soil.SW_mat_2, soil.SW_vec_2);

    return nothing
end




"""
    fit_soil_mat!(
                soil::SoilOpticals{FT},
                wls::WaveLengths{FT},
                ref_PAR::FT,
                ref_NIR::FT,
                method::TwoBandsFittingPoint
    ) where {FT<:AbstractFloat}

Fit soil albedo bands parameters, given
- `soil` [`SoilOpticals`](@ref) type structure
- `wls` [`WaveLengths`](@ref) type structure
- `ref_PAR` Target albedo for PAR
- `ref_NIR` Target albedo for NIR
- `method` [`TwoBandsFittingPoint`](@ref) type fitting method
"""
function fit_soil_mat!(
            soil::SoilOpticals{FT},
            wls::WaveLengths{FT},
            ref_PAR::FT,
            ref_NIR::FT,
            method::TwoBandsFittingPoint
) where {FT<:AbstractFloat}
    @unpack dry_NIR, dry_PAR, wet_NIR, wet_PAR = soil;

    # update soil PAR and NIR albedo
    soil.ρ_PAR = ref_PAR;
    soil.ρ_NIR = ref_NIR;

    # solve for weights
    soil.SW_vec_2[1] = (ref_PAR * wet_NIR - ref_NIR * wet_PAR) /
                       (wet_NIR * dry_PAR - wet_PAR * dry_NIR);
    soil.SW_vec_2[2] = (ref_PAR * dry_NIR - ref_NIR * dry_PAR) /
                       (dry_NIR * wet_PAR - dry_PAR * wet_NIR);

    # update vectors in soil
    mul!(soil.ρ_SW, soil.SW_mat_2, soil.SW_vec_2);

    return nothing
end
