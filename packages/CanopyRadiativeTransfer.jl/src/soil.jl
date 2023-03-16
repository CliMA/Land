# Soil albedo values from CLM
const SOIL_ALBEDOS = [0.36 0.61 0.25 0.50;
                      0.34 0.57 0.23 0.46;
                      0.32 0.53 0.21 0.42;
                      0.31 0.51 0.20 0.40;
                      0.30 0.49 0.19 0.38;
                      0.29 0.48 0.18 0.36;
                      0.28 0.45 0.17 0.34;
                      0.27 0.43 0.16 0.32;
                      0.26 0.41 0.15 0.30;
                      0.25 0.39 0.14 0.28;
                      0.24 0.37 0.13 0.26;
                      0.23 0.35 0.12 0.24;
                      0.22 0.33 0.11 0.22;
                      0.20 0.31 0.10 0.20;
                      0.18 0.29 0.09 0.18;
                      0.16 0.27 0.08 0.16;
                      0.14 0.25 0.07 0.14;
                      0.12 0.23 0.06 0.12;
                      0.10 0.21 0.05 0.10;
                      0.08 0.16 0.04 0.08];


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-09: rename function to soil_albedo!
#     2022-Jun-14: migrate the function from CanopyLayers
#     2022-Jun-14: add method to update broadband or hyperspectral soil albedo
#     2022-Jul-27: use albedo.α_CLM from ClimaCache v1.1.1, and remove option clm
#     2022-Jul-27: add albedo._θ control to HyperspectralSoilAlbedo method (fitting required)
#
#######################################################################################################################################################################################################
"""

    soil_albedo!(can::HyperspectralMLCanopy{FT}, soil::Soil{FT}) where {FT<:AbstractFloat}

Updates lower soil boundary reflectance, given
- `can` `HyperspectralMLCanopy` type struct
- `soil` `Soil` type struct

"""
function soil_albedo! end

soil_albedo!(can::HyperspectralMLCanopy{FT}, soil::Soil{FT}) where {FT<:AbstractFloat} = soil_albedo!(can, soil, soil.ALBEDO);

soil_albedo!(can::HyperspectralMLCanopy{FT}, soil::Soil{FT}, albedo::BroadbandSoilAlbedo{FT}) where {FT<:AbstractFloat} = (
    (; COLOR, LAYERS) = soil;
    @assert 1 <= COLOR <=20;

    # use CLM method
    if albedo.α_CLM
        _delta = max(0, FT(0.11) - FT(0.4) * LAYERS[1].θ);
        albedo.ρ_sw[1] = max(SOIL_ALBEDOS[COLOR,1], SOIL_ALBEDOS[COLOR,3] + _delta);
        albedo.ρ_sw[2] = max(SOIL_ALBEDOS[COLOR,2], SOIL_ALBEDOS[COLOR,4] + _delta);

        return nothing
    end;

    # use Yujie's method via relative soil water content
    _rwc = LAYERS[1].θ / LAYERS[1].VC.Θ_SAT;
    albedo.ρ_sw[1] = SOIL_ALBEDOS[COLOR,1] * (1 - _rwc) + _rwc * SOIL_ALBEDOS[COLOR,3];
    albedo.ρ_sw[2] = SOIL_ALBEDOS[COLOR,2] * (1 - _rwc) + _rwc * SOIL_ALBEDOS[COLOR,4];

    return nothing
);

soil_albedo!(can::HyperspectralMLCanopy{FT}, soil::Soil{FT}, albedo::HyperspectralSoilAlbedo{FT}) where {FT<:AbstractFloat} = (
    (; WLSET) = can;
    (; COLOR, LAYERS) = soil;
    @assert 1 <= COLOR <=20;

    # if the change of swc is lower than 0.01, do nothing
    if abs(LAYERS[1].θ - albedo._θ) < 0.01
        return nothing
    end;

    # use CLM method or Yujie's method
    _rwc::FT = LAYERS[1].θ / LAYERS[1].VC.Θ_SAT;
    _par::FT = SOIL_ALBEDOS[COLOR,1] * (1 - _rwc) + _rwc * SOIL_ALBEDOS[COLOR,3];
    _nir::FT = SOIL_ALBEDOS[COLOR,2] * (1 - _rwc) + _rwc * SOIL_ALBEDOS[COLOR,4];

    if albedo.α_CLM
        _delta = max(0, FT(0.11) - FT(0.4) * LAYERS[1].θ);
        _par = max(SOIL_ALBEDOS[COLOR,1], SOIL_ALBEDOS[COLOR,3] + _delta);
        _nir = max(SOIL_ALBEDOS[COLOR,2], SOIL_ALBEDOS[COLOR,4] + _delta);
    end;

    # if fitting is disabled, use broadband directly
    if !albedo.FITTING
        albedo.ρ_sw[WLSET.IΛ_PAR] .= _par;
        albedo.ρ_sw[WLSET.IΛ_NIR] .= _nir;

        return nothing
    end;

    # make an initial guess of the weights
    albedo._ρ_sw[WLSET.IΛ_PAR] .= _par;
    albedo._ρ_sw[WLSET.IΛ_NIR] .= _nir;
    albedo._weight .= pinv(albedo.MAT_ρ) * albedo._ρ_sw;

    # function to solve for weights
    @inline _fit(x::Vector{FT}) where {FT<:AbstractFloat} = (
        mul!(albedo._ρ_sw, albedo.MAT_ρ, x);
        albedo._tmp_vec_nir .= abs.(view(albedo._ρ_sw,WLSET.IΛ_NIR) .- _nir);
        _diff = ( mean( view(albedo._ρ_sw,WLSET.IΛ_PAR) ) - _par ) ^ 2 + mean( albedo._tmp_vec_nir ) ^ 2;

        return -_diff
    );

    # solve for weights
    _ms  = ReduceStepMethodND{FT}(x_mins = FT[-2,-2,-2,-2], x_maxs = FT[2,2,2,2], x_inis = albedo._weight, Δ_inis = FT[0.1,0.1,0.1,0.1]);
    _tol = SolutionToleranceND{FT}(FT[0.001,0.001,0.001,0.001], 50);
    _sol = find_peak(_fit, _ms, _tol);
    albedo._weight .= _sol;

    # update vectors in soil
    mul!(albedo.ρ_sw, albedo.MAT_ρ, albedo._weight);

    # update the albedo._θ to avoid calling this function too many times
    albedo._θ = LAYERS[1].θ;

    return nothing
);
