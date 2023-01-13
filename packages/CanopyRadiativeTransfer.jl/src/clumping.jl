#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-02: migrate the function from CanopyLayers
#     2022-Jun-02: rename the function from clumping_factor! to clumping_index!
#
#######################################################################################################################################################################################################
"""

    clumping_index!(can::HyperspectralMLCanopy, angles::SunSensorGeometry{FT}) where {FT<:AbstractFloat}

Update the clumping index, given
- `can` `HyperspectralMLCanopy` type canopy
- `angles` `SunSensorGeometry` type angles

"""
function clumping_index!(can::HyperspectralMLCanopy, angles::SunSensorGeometry{FT}) where {FT<:AbstractFloat}
    (; Ω_A, Ω_B) = can;

    can.ci = Ω_A + Ω_B * (1 - cosd(angles.sza));

    return nothing
end
