###############################################################################
#
# Update canopy matrices
#
###############################################################################
"""
    canopy_matrices!(
                leaves::Array{LeafBios{FT},1},
                can_opt::CanopyOpticals{FT}
    ) where {FT<:AbstractFloat}

Compute scattering coefficient matrices for direct and diffuse light given
    geometry dependent overall extinction coefficients and pigment dependent
    leaf reflectance and transmission (computed via fluspect). This function
    has to be called before [`short_wave!`](@ref) can be used.
- `leaves` Array of [`LeafBios`](@ref) type struct
- `can_opt` [`CanopyOpticals`](@ref) type struct
"""
function canopy_matrices!(
            leaves::Array{LeafBios{FT},1},
            can_opt::CanopyOpticals{FT}
) where {FT<:AbstractFloat}
    # 1. unpack values
    @unpack ddb, ddf, dob, dof, sdb, sdf, sob, sof = can_opt;

    # 2. Calculation of reflectance
    nLayer = size(can_opt.sigb)[2];
    @inbounds for i=1:nLayer
        if length(leaves)>1
            τ_SW = leaves[i].τ_SW;
            ρ_SW = leaves[i].ρ_SW;
        else
            τ_SW = leaves[1].τ_SW;
            ρ_SW = leaves[1].ρ_SW;
        end

        #CF: Right now, canopy geometry is the same everywhere,
        #    can be easily extended to layers as well.
        @inbounds for j=1:size(can_opt.sigb, 1)
            # diffuse backscatter coefficient for diffuse incidence
            can_opt.sigb[j,i] = ddb * ρ_SW[j] + ddf * τ_SW[j];
            # diffuse forwardscatter coefficient for diffuse incidence
            can_opt.sigf[j,i] = ddf * ρ_SW[j] + ddb * τ_SW[j];
            # diffuse backscatter coefficient for specular incidence
            can_opt.sb[j,i]   = sdb * ρ_SW[j] + sdf * τ_SW[j];
            # diffuse forwardscatter coefficient for specular incidence
            can_opt.sf[j,i]   = sdf * ρ_SW[j] + sdb * τ_SW[j];
            # directional backscatter  coefficient for diffuse incidence
            can_opt.vb[j,i]   = dob * ρ_SW[j] + dof * τ_SW[j];
            # directional forwardscatter coefficient for diffuse incidence
            can_opt.vf[j,i]   = dof * ρ_SW[j] + dob * τ_SW[j];
            # bidirectional scattering coefficent (directional-directional)
            can_opt.w[j,i]    = sob * ρ_SW[j] + sof * τ_SW[j];
        end
    end

    # 3. attenuation
    can_opt.a .= 1 .- can_opt.sigf;

    return nothing
end
