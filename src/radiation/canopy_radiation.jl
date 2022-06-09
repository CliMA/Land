#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-09: migrate CanopyRads as CanopyRadiationProfile
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to store canopy radiation profiles

# Fields

$(TYPEDFIELDS)

"""
mutable struct CanopyRadiationProfile{FT<:AbstractFloat}
    # diagnostic variables that change with time
    "Downwelling diffuse short-wave radiation at each canopy layer boundary `[mW m⁻² nm⁻¹]`"
    e_diffuse_down::Matrix{FT}
    "Upwelling diffuse short-wave radiation at each canopy layer boundary `[mW m⁻² nm⁻¹]`"
    e_diffuse_up::Matrix{FT}
    "Solar directly radiation at each canopy layer boundary `[mW m⁻² nm⁻¹]`"
    e_direct::Matrix{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-09: add constructor
#
#######################################################################################################################################################################################################
"""

    CanopyRadiationProfile{FT}(; n_layer::Int = 20, n_λ::Int = 114) where {FT<:AbstractFloat}

Construct a struct to store canopy radiation profiles, given
- `n_layer` Number of canopy layers
- `n_λ` Number of wavelength bins
"""
CanopyRadiationProfile{FT}(; n_layer::Int = 20, n_λ::Int = 114) where {FT<:AbstractFloat} = (
    return CanopyOpticalProperty{FT}(
                zeros(FT,n_λ,n_layer+1),    # e_diffuse_down
                zeros(FT,n_λ,n_layer+1),    # e_diffuse_up
                zeros(FT,n_λ,n_layer+1)     # e_direct
    )
);
