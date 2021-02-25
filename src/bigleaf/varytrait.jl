###############################################################################
#
# Vary the traits in the SPACSimple struct
#
###############################################################################
"""
    vary_spac!(node::SPACSimple{FT},
               weat_years::DataFrame,
               factor_to_vary::String,
               ratio_to_vary::FT
    ) where {FT<:AbstractFloat}

Generalized function to vary the SoilPlantAirContinuum to run sensitivity
    analysis, given
- `node` SPACSimple type struct
- `weat_years` Weather data
- `factor_to_vary` Which parameter to vary
- `ratio_to_vary` How much to vary
"""
function vary_spac!(
            node::SPACSimple{FT},
            weat_years::DataFrame,
            factor_to_vary::String,
            ratio_to_vary::FT
) where {FT<:AbstractFloat}
    # choose which parameter(s) to vary
    if factor_to_vary=="wb"
        (node.hs.root.vc).b *= ratio_to_vary;
        (node.hs.stem.vc).b *= ratio_to_vary;
        (node.hs.leaf.vc).b *= ratio_to_vary;
    elseif factor_to_vary=="wc"
        (node.hs.root.vc).c *= ratio_to_vary;
        (node.hs.stem.vc).c *= ratio_to_vary;
        (node.hs.leaf.vc).c *= ratio_to_vary;
    elseif factor_to_vary=="wk"
        node.hs.root.k_max *= ratio_to_vary;
        node.hs.stem.k_max *= ratio_to_vary;
        node.hs.leaf.k_sla *= ratio_to_vary;
        node.hs.root.k_element .*= ratio_to_vary;
        node.hs.stem.k_element .*= ratio_to_vary;
        node.hs.leaf.k_element .*= ratio_to_vary;
    elseif factor_to_vary=="kw"
        node.hs.root.k_max *= ratio_to_vary;
        node.hs.stem.k_max *= ratio_to_vary;
        node.hs.root.k_element .*= ratio_to_vary;
        node.hs.stem.k_element .*= ratio_to_vary;
    elseif factor_to_vary=="kl"
        node.hs.leaf.k_sla *= ratio_to_vary;
        node.hs.leaf.k_element .*= ratio_to_vary;
    elseif factor_to_vary=="cv"
        node.c_vmax *= ratio_to_vary;
    elseif factor_to_vary=="cc"
        node.c_cons *= ratio_to_vary;
    elseif factor_to_vary=="gm"
        node.g_max *= ratio_to_vary;
    elseif factor_to_vary=="ga"
        node.gaba *= ratio_to_vary;
        node.lai  /= ratio_to_vary;
    elseif factor_to_vary=="sd"
        node.h_soil *= ratio_to_vary;
        node.hs.root.p_gravity .*= ratio_to_vary;
    elseif factor_to_vary=="ca"
        (weat_years).CO2 *= ratio_to_vary;
    elseif factor_to_vary=="rh"
        (weat_years).D *= ratio_to_vary;
    elseif factor_to_vary=="ta"
        _multip1 = saturation_vapor_pressure.((weat_years).Tair .+ FT(273.15));
        (weat_years).Tair .+= ratio_to_vary;
        _multip2 = saturation_vapor_pressure.((weat_years).Tair .+ FT(273.15));
        (weat_years).D    .*= _multip2 ./ _multip1;
    else
        @warn "Invalid parameter provided, nothing has been changed";
    end

    return nothing
end
