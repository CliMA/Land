#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-30: migrate function from older version StomataModels.jl
#     2022-Jul-01: add method to tune stomatal opening based on relative hydraulic conductance at leaf xylem end
#     2022-Jul-01: add method to tune stomatal opening based on relative hydraulic conductance of the soil
#     2022-Jul-01: add method to tune stomatal opening based on soil potential or leaf pressure
#     2022-Jul-01: fix a typo in function call
#     2022-Jul-11: deflate documentations
#     2022-Jul-12: move function from StomataModels.jl to PlantHydraulics.jl
#
#######################################################################################################################################################################################################
"""

    β_factor(f::Function, vc::AbstractXylemVC{FT}, x_25::FT) where {FT<:AbstractFloat}
    β_factor(f::Function, vc::AbstractSoilVC{FT}, x_25::FT) where {FT<:AbstractFloat}
    β_factor(f::Function, x_25::FT) where {FT<:AbstractFloat}

Return the β factor based on relative conductance or soil potential/pressure, given
- `f` Function to translate relative k to β, for example f(x) = x, f(x) = x², and f(x) = sqrt(x) for x in [0,1]
- `vc` Leaf vulnerability curve or soil vulnerability curve (moisture retention curve)
- `x_25` Leaf xylem pressure corrected to 25 °C, soil water potential corrected to 25 °C (forcing on roots, note that this function may not be useful for plants with salt stress), or soil water
    content

"""
function β_factor end

β_factor(f::Function, vc::AbstractXylemVC{FT}, x_25::FT) where {FT<:AbstractFloat} = FT(f(relative_hydraulic_conductance(vc, x_25)));

β_factor(f::Function, vc::AbstractSoilVC{FT}, x_25::FT) where {FT<:AbstractFloat} = FT(f(relative_hydraulic_conductance(vc, true, x_25)));

β_factor(f::Function, x_25::FT) where {FT<:AbstractFloat} = FT(f(x_25));


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-12: add function to update beta factor for empirical models
#     2022-Jul-12: add methods for MonoElementSPAC
#
#######################################################################################################################################################################################################
"""
This function updates the beta factor for SPAC if empirical models are used. The method is meant to support all SPAC defined in ClimaCache.jl:
- `MonoElementSPAC`
- `MonoMLGrassSPAC`
- `MonoMLPalmSPAC`
- `MonoMLTreeSPAC`

"""
function β_factor! end

"""

    β_factor!(spac::MonoElementSPAC{FT}) where {FT<:AbstractFloat}

Update the beta factor for the LEAF component in SPAC, given
- `spac` `MonoElementSPAC` type SPAC

"""
β_factor!(spac::MonoElementSPAC{FT}) where {FT<:AbstractFloat} = β_factor!(spac, spac.LEAF.SM);

β_factor!(spac::MonoElementSPAC{FT}, sm::Union{AndereggSM{FT}, EllerSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}}) where {FT<:AbstractFloat} = nothing;

β_factor!(spac::MonoElementSPAC{FT}, sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}) where {FT<:AbstractFloat} = β_factor!(spac, sm.Β);

β_factor!(spac::MonoElementSPAC{FT}, β::BetaFunction{FT}) where {FT<:AbstractFloat} = β_factor!(spac, β, β.PARAM_X);

β_factor!(spac::MonoElementSPAC{FT}, β::BetaFunction{FT}, param_x::BetaParameterKleaf) where {FT<:AbstractFloat} = (
    _f_st = relative_surface_tension(spac.LEAF.t);

    β.β₁ = β_factor(β.FUNC, spac.LEAF.HS.VC, spac.LEAF.HS.p_element[end] / _f_st);

    return nothing
);

β_factor!(spac::MonoElementSPAC{FT}, β::BetaFunction{FT}, param_x::BetaParameterKsoil) where {FT<:AbstractFloat} = (
    _f_st = relative_surface_tension(spac.ROOT.t);

    β.β₁ = β_factor(β.FUNC, spac.ROOT.HS.SH, spac.ROOT.HS.p_ups / _f_st);

    return nothing
);

β_factor!(spac::MonoElementSPAC{FT}, β::BetaFunction{FT}, param_x::BetaParameterPleaf) where {FT<:AbstractFloat} = (
    _f_st = relative_surface_tension(spac.LEAF.t);

    β.β₁ = β_factor(β.FUNC, spac.LEAF.HS.p_element[end] / _f_st);

    return nothing
);

β_factor!(spac::MonoElementSPAC{FT}, β::BetaFunction{FT}, param_x::BetaParameterPsoil) where {FT<:AbstractFloat} = (
    _f_st = relative_surface_tension(spac.ROOT.t);

    β.β₁ = β_factor(β.FUNC, spac.ROOT.HS.p_ups / _f_st);

    return nothing
);

β_factor!(spac::MonoElementSPAC{FT}, β::BetaFunction{FT}, param_x::BetaParameterΘ) where {FT<:AbstractFloat} = (
    _f_st = relative_surface_tension(spac.ROOT.t);

    β.β₁ = β_factor(β.FUNC, soil_θ(spac.ROOT.HS.SH, spac.ROOT.HS.p_ups / _f_st));

    return nothing
);
