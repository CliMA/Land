###############################################################################
#
# Pre-defined parameter sets
#
###############################################################################
"""
    create_VanGenuchten(
                name::String,
                α::FT,
                n::FT,
                Θs::FT,
                Θr::FT) where {FT<:AbstractFloat}

Create a [`AbstractSoilVC`](@ref) type of soil VC, given
- `name` Soil type name
- `α` Soil α
- `n` Soil n
- `Θs` SWC at saturation
- `Θr` Residual SWC
"""
function create_VanGenuchten(
            name::String,
            α::FT,
            n::FT,
            Θs::FT,
            Θr::FT
) where {FT<:AbstractFloat}
    return VanGenuchten{FT}(stype = name,
                            α     = α   ,
                            n     = n   ,
                            Θs    = Θs  ,
                            Θr    = Θr  )
end




"""
Create VanGenuchten VC for sand soil
"""
VGSand(FT) = create_VanGenuchten(
                        "Sand",
                        FT(1479.5945),
                        FT(2.68),
                        FT(0.43),
                        FT(0.045))

"""
Create VanGenuchten VC for loamy sand soil
"""
VGLoamySand(FT) = create_VanGenuchten(
                        "Loamy Sand",
                        FT(1265.3084),
                        FT(2.28),
                        FT(0.41),
                        FT(0.057))

"""
Create VanGenuchten VC for Sandy Loam soil
"""
VGSandyLoam(FT) = create_VanGenuchten(
                        "Sandy Loam",
                        FT(765.3075),
                        FT(1.89),
                        FT(0.41),
                        FT(0.065))

"""
Create VanGenuchten VC for Loamy Sand 2 soil
"""
VGLoamySand2(FT) = create_VanGenuchten(
                        "Loamy Sand 2",
                        FT(367.3476),
                        FT(1.56),
                        FT(0.43),
                        FT(0.078))

"""
Create VanGenuchten VC for Sandy Clay Loam soil
"""
VGSandyClayLoam(FT) = create_VanGenuchten(
                        "Sandy Clay Loam",
                        FT(602.0419),
                        FT(1.48),
                        FT(0.39),
                        FT(0.1))

"""
Create VanGenuchten VC for Silty Loam soil
"""
VGSiltyLoam(FT) = create_VanGenuchten(
                        "Silty Loam",
                        FT(204.082),
                        FT(1.41),
                        FT(0.45),
                        FT(0.067))

"""
Create VanGenuchten VC for Silt soil
"""
VGSilt(FT) = create_VanGenuchten(
                        "Silt",
                        FT(163.2656),
                        FT(1.37),
                        FT(0.46),
                        FT(0.034))

"""
Create VanGenuchten VC for Clay Loam soil
"""
VGClayLoam(FT) = create_VanGenuchten(
                        "Clay Loam",
                        FT(193.8779),
                        FT(1.31),
                        FT(0.41),
                        FT(0.095))

"""
Create VanGenuchten VC for Silty Clay Loam soil
"""
VGSiltyClayLoam(FT) = create_VanGenuchten(
                        "Silty Clay Loam",
                        FT(102.041),
                        FT(1.23),
                        FT(0.43),
                        FT(0.089))

"""
Create VanGenuchten VC for Sandy Clay Loam 2 soil
"""
VGSandyClayLoam2(FT) = create_VanGenuchten(
                        "Sandy Clay Loam 2",
                        FT(275.5107),
                        FT(1.23),
                        FT(0.38),
                        FT(0.1))

"""
Create VanGenuchten VC for Silty Clay soil
"""
VGSiltyClay(FT) = create_VanGenuchten(
                        "Silty Clay",
                        FT(51.0205),
                        FT(1.09),
                        FT(0.36),
                        FT(0.07))

"""
Create VanGenuchten VC for Clay soil
"""
VGClay(FT) = create_VanGenuchten(
                        "Clay",
                        FT(81.6328),
                        FT(1.09),
                        FT(0.38),
                        FT(0.068))
