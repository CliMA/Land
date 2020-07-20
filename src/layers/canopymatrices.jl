###############################################################################
#
# Update canopy matrices
#
###############################################################################
"""
    canopy_matrices!(leaf_array::Array{LeafBios{FT},1}, can_opt::CanopyOpticals{FT}) where {FT<:AbstractFloat}

Compute scattering coefficient matrices for direct and diffuse light given
    geometry dependent overall extinction coefficients and pigment dependent
    leaf reflectance and transmission (computed via fluspect). This function
    has to be called before `simulate_short_wave!` can be used.
- `leaf_array` An array of [`LeafBios`](@ref) type struct (i.e. leaf optical
    properties can change with canopy height)
- `can_opt` An [`CanopyOpticals`](@ref) type struct, will be updated within
    this function call.
"""
function canopy_matrices!(
            leaf_array::Array{LeafBios{FT},1},
            can_opt::CanopyOpticals{FT}
) where {FT<:AbstractFloat}
    # 2. Calculation of reflectance
    # 2.1 reflectance, transmittance factors in a thin layer the following are vectors with length [nl,nwl]
    nlayers = size(can_opt.sigb)[2];
    @inbounds for i=1:nlayers
        if length(leaf_array)>1
            τ_SW = leaf_array[i].τ_SW;
            ρ_SW = leaf_array[i].ρ_SW;
        else
            τ_SW = leaf_array[1].τ_SW;
            ρ_SW = leaf_array[1].ρ_SW;
        end
        #CF: Right now, canopy geometry is the same everywhere, can be easily extended to layers as well.
        @inbounds for j=1:size(can_opt.sigb, 1)
            # [nl,nwl] diffuse backscatter scattering coefficient for diffuse incidence
            can_opt.sigb[j,i] = can_opt.ddb * ρ_SW[j] + can_opt.ddf * τ_SW[j];
            # [nl,nwl] diffuse forward scattering coefficient for diffuse incidence
            can_opt.sigf[j,i] = can_opt.ddf * ρ_SW[j] + can_opt.ddb * τ_SW[j];
            # [nl,nwl] diffuse backscatter scattering coefficient for specular incidence
            can_opt.sb[j,i]   = can_opt.sdb * ρ_SW[j] + can_opt.sdf * τ_SW[j];
            # [nl,nwl] diffuse forward scattering coefficient for specular incidence
            can_opt.sf[j,i]   = can_opt.sdf * ρ_SW[j] + can_opt.sdb * τ_SW[j];
            # [nl,nwl] directional backscatter scattering coefficient for diffuse incidence
            can_opt.vb[j,i]   = can_opt.dob * ρ_SW[j] + can_opt.dof * τ_SW[j];
            # [nl,nwl] directional forward scattering coefficient for diffuse incidence
            can_opt.vf[j,i]   = can_opt.dof * ρ_SW[j] + can_opt.dob * τ_SW[j];
            # [nl,nwl] bidirectional scattering coefficent (directional-directional)
            can_opt.w[j,i]    = can_opt.sob * ρ_SW[j] + can_opt.sof * τ_SW[j];
        end
    end

    # [nl, nwl] attenuation
    can_opt.a .= 1 .- can_opt.sigf;

    return nothing
end
