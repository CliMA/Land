"""
    update_canopy_from_rt_module!(tree::Tree, canopy_rt::struct_canopy, canOpt_rt::struct_canopyOptProps, canRad_rt::struct_canopyRadiation)

Updates canopy information from the RT module, given
- `tree` A [`Tree`](@ref) type
- `canopy_rt` A [`struct_canopy`] type from the Plant module
- `cadOpt_rt` A [`struct_canopyOptProps`] type from the CanopyRT module
- `canRad_rt` A [`struct_canopyRadiation`] type from the CanopyRT module

This interface function is pending for leaf temperature...
"""
function update_canopy_from_rt_module!(tree::Tree, canopy_rt::struct_canopy, canOpt_rt::struct_canopyOptProps, canRad_rt::struct_canopyRadiation)
    # fraction of sunlit leaves in each layer
    fraction_sl = repeat(canopy_rt.lidf, outer=[ length(canopy_rt.lazitab) ]) / length(canopy_rt.lazitab)

    # update the PAR from canopyRT module to the Plant module
    if tree.canopy.n_layer == canopy_rt.nlayers
        nlayers = canopy_rt.nlayers
        for pl_layer in 1:nlayers
            canopyi  = tree.canopy.canopy_list[pl_layer]
            rt_layer = nlayers + 1 - pl_layer
            
            # set the diffuse par to all the leaves
            canopyi.par_list .= canRad_rt.absPAR_shadeCab[rt_layer] * 1E6
            # add the direct PAR to sunlit leaves
            canopyi.par_list[1:end-1] .+= reshape(canRad_rt.absPAR_sunCab[:,:,rt_layer],(:,1))[:,1] * 1E6

            # calculate the fraction of sunlit and shaded leaves
            f_view = (canOpt_rt.Ps[rt_layer]+canOpt_rt.Ps[rt_layer+1]) / 2
            la_new = canopyi.la .* [f_view .* fraction_sl; 1-f_view]
            canopyi.f_view  = f_view
            canopyi.la_list = la_new

            # update leaf temperature
            canopyi.t_list[1:end-1] .= reshape(canRad_rt.T_sun3D[:,:,rt_layer],(:,1))
            canopyi.t_list[end]      = canRad_rt.T_shade[rt_layer]
        end
    else
        println("Error: the canopy layer in the Tree differs from that in the struct_canopy!")
    end
end
