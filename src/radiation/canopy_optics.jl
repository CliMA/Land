#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-07: add CanopyOptics struct (will be a field for canopy structure)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for Verhoef LIDF algorithm

# Fields

$(TYPEDFIELDS)

"""
mutable struct CanopyOpticalProperty{FT<:AbstractFloat}
    # diagnostic variables that change with time
    "Diffuse -> Diffuse backscatter weight"
    ddb::FT
    "Diffuse -> Diffuse forward scatter weight"
    ddf::FT
    "Diffuse -> Outgoing backscatter weight"
    dob::FT
    "Diffuse -> Outgoing forward scatter weight"
    dof::FT
    "Outgoing beam extinction coefficient weight"
    ko ::FT
    "Solar beam extinction coefficient weight"
    ks ::FT
    "Solar -> Diffuse backscatter weight"
    sdb::FT
    "Solar -> Diffuse forward scatter weight"
    sdf::FT
    "Solar -> Outgoing weight of specular2directional backscatter coefficient"
    sob::FT
    "Solar -> Outgoing weight of specular2directional forward coefficient"
    sof::FT

    # caches to speed up calculations
    "Weighted sum of cos²(inclination)"
    _bf::FT
    "Outgoing beam extinction coefficient weights at different inclination angles"
    _ko::Vector{FT}
    "Solar beam extinction coefficient weights at different inclination angles"
    _ks::Vector{FT}
    "Backward scattering coefficients at different inclination angles"
    _sb::Vector{FT}
    "Forward scattering coefficients at different inclination angles"
    _sf::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-07: add constructor
#
#######################################################################################################################################################################################################
"""

    CanopyOpticalProperty{FT}(; n_incl::Int = 9) where {FT<:AbstractFloat}

Construct a struct to store canopy optical properties
"""
CanopyOpticalProperty{FT}(; n_incl::Int = 9) where {FT<:AbstractFloat} = (
    return CanopyOpticalProperty{FT}(
                0,                  # ddb
                0,                  # ddf
                0,                  # dob
                0,                  # dof
                0,                  # ko
                0,                  # ks
                0,                  # sdb
                0,                  # sdf
                0,                  # sob
                0,                  # sof
                0,                  # _bf
                zeros(FT,n_incl),   # _ko
                zeros(FT,n_incl),   # _ks
                zeros(FT,n_incl),   # _sb
                zeros(FT,n_incl)    # _sf
    )
);
