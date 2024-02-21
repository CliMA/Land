###############################################################################
#
# Soil optical parameters
#
###############################################################################
"""
    mutable struct SoilOpticals{FT}

A struct of soil optical parameters

# Fields
$(TYPEDFIELDS)
"""
mutable struct SoilOpticals{FT}
    "Soil surface temperature"
    T::FT
    "Soil color class"
    color::Int
    "Whether to expand from broadband to hyperspectral"
    hyperspectral::Bool

    # broad band showrtwave albedo
    "Shortwave albedo for NIR"
    ρ_NIR::FT
    "Shortwave albedo for PAR"
    ρ_PAR::FT

    # hyperspectral shortwave albedo
    "Shortwave albedo that matches `WL` from [`WaveLengths`](@ref)"
    ρ_SW::Vector{FT}
    "Shortwave albedo that matches `WLF` from [`WaveLengths`](@ref)"
    ρ_SW_SIF::Vector{FT}
    "Shortwave absorption that equals `1 - ρ_SW`"
    ε_SW::Vector{FT}

    # hyperspectral soil albedo
    "Shortwave albedo matrix from 4 bands with `WL` from [`WaveLengths`](@ref)"
    SW_mat_4::Matrix{FT}
    "Shortwave albedo matrix from 2 bands with `WL` from [`WaveLengths`](@ref)"
    SW_mat_2::Matrix{FT}
    "Shortwave albedo weight from 4 bands"
    SW_vec_4::Vector{FT}
    "Shortwave albedo weight from 2 bands"
    SW_vec_2::Vector{FT}

    # hyperspectral longwave albedo
    "Longtwave albedo"
    ρ_LW::Vector{FT}

    # cache that stores mean band values (used to speed up calculations)
    "Mean value for day band 1 in NIR region"
    dry_NIR::FT
    "Mean value for day band 1 in PAR region"
    dry_PAR::FT
    "Mean value for day band 1 in NIR region"
    wet_NIR::FT
    "Mean value for day band 1 in PAR region"
    wet_PAR::FT
end


# constructor
SoilOpticals(wls::WaveLengths{FT}) where {FT<:AbstractFloat} = (
    (; nWL, nWLF, opti_file) = wls;

    # rescale the data to match the steppings of wavelength set
    _res_4 = FT[read_nc(opti_file, "GSV_1") read_nc(opti_file, "GSV_2") read_nc(opti_file, "GSV_3") read_nc(opti_file, "GSV_4")];
    _res_2 = FT[_res_4[:,1] _res_4[:,end]];

    return SoilOpticals{FT}(T₂₅(FT), 1, true, FT(0.2), FT(0.2), ones(FT,nWL).*FT(0.2),
                            ones(FT,nWLF)*FT(0.2), ones(FT,nWL).*FT(0.8),
                            _res_4, _res_2, ones(FT,4), ones(FT,2),
                            FT[0.1], FT(0.5), FT(0.5), FT(0.5), FT(0.5))
);
