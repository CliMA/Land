###############################################################################
#
# Big leaf canopy model partition
#
###############################################################################
"""
    big_leaf_partition!(
                node::SPACSimple{FT},
                zenith::FT,
                r_all::FT
    ) where {FT <:AbstractFloat}

Partition the big-leaf canopy into sunlit and shaded layers, given
- `partition` Container for partition
- `zenith` Zenith angle in degree
- `r_all` Total radiation in `[W m⁻²]`
"""
function big_leaf_partition!(
            node::SPACSimple{FT},
            zenith::FT,
            r_all::FT
) where {FT <:AbstractFloat}
    # unpack values
    @unpack laba, lai = node;

    # 1. use big_leaf_partition function from CanopyLayers module
    ratio, q_slm, q_shm, e_sl, e_sh = big_leaf_partition(lai, zenith, r_all);

    # 2. update the information
    node.container2L.frac_sl = ratio;
    node.container2L.frac_sh = 1 - ratio;
    node.container2L.la_sl   = laba * ratio;
    node.container2L.la_sh   = laba * (1 - ratio);
    node.container2L.lai_sl  = lai * ratio;
    node.container2L.lai_sh  = lai * (1 - ratio);
    node.container2L.par_sl  = q_slm;
    node.container2L.par_sh  = q_shm;
    node.container2L.rad_sl  = e_sl;
    node.container2L.rad_sh  = e_sh;

    return nothing
end
