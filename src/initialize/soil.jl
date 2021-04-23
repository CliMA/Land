###############################################################################
#
# Pre-defined parameter sets
#
###############################################################################
"""
    create_soil_VC(
                vc::AbstractSoilVC{FT},
                name::String,
                α::Number,
                n::Number,
                Θs::Number,
                Θr::Number
    ) where {FT<:AbstractFloat}
    create_soil_VC(
                vc::AbstractSoilVC{FT},
                name::String
    ) where {FT<:AbstractFloat}

Create a [`AbstractSoilVC`](@ref) type of soil VC, given
- `vc` [`AbstractSoilVC`](@ref) type identifier
- `name` Soil type name. Supported names include
  - "Sand"
  - "Loamy Sand"
  - "Sandy Loam"
  - "Loam"
  - "Sandy Clay Loam"
  - "Silt Loam"
  - "Silt"
  - "Clay Loam"
  - "Silty Clay Loam"
  - "Sandy Clay"
  - "Silty Clay"
  - "Clay"
- `α` Soil α
- `n` Soil n
- `Θs` SWC at saturation
- `Θr` Residual SWC
"""
function create_soil_VC(
            vc::BrooksCorey{FT},
            name::String,
            α::Number,
            n::Number,
            Θs::Number,
            Θr::Number
) where {FT<:AbstractFloat}
    _vc_vG = VanGenuchten{FT}(stype = name,
                                  α = α,
                                  n = n,
                                 Θs = Θs,
                                 Θr = Θr)
    _vc_BC = BrooksCorey{FT}(stype = name,
                                 b = (2n-1)/(n-1),
                                ϕs = 1/α,
                                Θs = Θs,
                                Θr = Θr)
    fit_soil_VC!(_vc_vG, _vc_BC);

    return _vc_BC
end




function create_soil_VC(
            vc::VanGenuchten{FT},
            name::String,
            α::Number,
            n::Number,
            Θs::Number,
            Θr::Number
) where {FT<:AbstractFloat}
    return VanGenuchten{FT}(stype = name,
                                α = α,
                                n = n,
                               Θs = Θs,
                               Θr = Θr)
end




function create_soil_VC(
            vc::AbstractSoilVC{FT},
            name::String
) where {FT<:AbstractFloat}
    # Parameters from Silt soil
    paras = [ 163.2656, 1.37, 0.46, 0.034];

    if name=="Sand"
        paras = [1479.5945, 2.68, 0.43, 0.045];
    elseif name=="Loamy Sand"
        paras = [1265.3084, 2.28, 0.41, 0.057];
    elseif name=="Sandy Loam"
        paras = [ 765.3075, 1.89, 0.41, 0.065];
    elseif name=="Loam"
        paras = [ 367.3476, 1.56, 0.43, 0.078];
    elseif name=="Sandy Clay Loam"
        paras = [ 602.0419, 1.48, 0.39, 0.100];
    elseif name=="Silt Loam"
        paras = [ 204.0820, 1.41, 0.45, 0.067];
    elseif name=="Silt"
        paras = [ 163.2656, 1.37, 0.46, 0.034];
    elseif name=="Clay Loam"
        paras = [ 193.8779, 1.31, 0.41, 0.095];
    elseif name=="Silty Clay Loam"
        paras = [ 102.0410, 1.23, 0.43, 0.089];
    elseif name== "Sandy Clay"
        paras = [ 275.5107, 1.23, 0.38, 0.100];
    elseif name=="Silty Clay"
        paras = [  51.0205, 1.09, 0.36, 0.070];
    elseif name=="Clay"
        paras = [  81.6328, 1.09, 0.38, 0.068];
    else
        @warn twarn("Soil type $(name) not recognized, use Silt instead.");
    end

    return create_soil_VC(vc, name, paras...)
end
