###############################################################################
#
# Update soil from water content or pressure
#
###############################################################################
"""
    soil_moisture_swc!(node::SPACSimple{FT}, swc::FT) where {FT<:AbstractFloat}

Update soil moisture and soil matrix potential, given
- `node` [`SPACSimple`](@ref) type struct
- `swc` Given soil water content
"""
function soil_moisture_swc!(
            node::SPACSimple{FT},
            swc::FT
) where {FT<:AbstractFloat}
    if swc < node.mswc
        node.swc = swc;
    else
        node.swc = node.mswc;
    end

    node.p_soil = soil_p_25_swc(node.hs.root.sh, node.swc) * node.hs.root.f_st;

    return nothing
end




"""
    soil_moisture_p!(node::SPACSimple{FT}, p::FT) where {FT<:AbstractFloat}

Update soil moisture and soil matrix potential, given
- `node` [`SPACSimple`](@ref) type struct
- `p` Given soil maxtrix potential
"""
function soil_moisture_p!(
            node::SPACSimple{FT},
            p::FT
) where {FT<:AbstractFloat}
    node.p_soil = p;
    p_25        = p / node.hs.root.f_st;
    node.swc    = soil_swc(node.hs.root.sh, p_25);

    return nothing
end




"""
    soil_moisture_p25!(
                node::SPACSimple{FT},
                p_25::FT
    ) where {FT<:AbstractFloat}

Update soil moisture and soil matrix potential, given
- `node` [`SPACSimple`](@ref) type struct
- `p_25` Given soil maxtrix potential at 25 Celcius
"""
function soil_moisture_p25!(
            node::SPACSimple{FT},
            p_25::FT
) where {FT<:AbstractFloat}
    node.p_soil = p_25 * node.hs.root.f_st;
    node.swc    = soil_swc(node.hs.root.sh, p_25);

    return nothing
end




"""
    soil_moisture!(
                node::SPACSimple{FT},
                flow::FT,
                Δt::FT = FT(1)
    ) where {FT<:AbstractFloat}

Update soil moisture and soil matrix potential, given
- `node` [`SPACSimple`](@ref) type struct
- `flow` Mean outlet flow rate in `[Kg h⁻¹]`
- `Δt` Time period in `[h]`
"""
function soil_moisture!(
            node::SPACSimple{FT},
            flow::FT,
            Δt::FT = FT(1)
) where {FT<:AbstractFloat}
    # 1. positive flow means out, flow in Kg h⁻¹, dt in h
    m_all  = node.swc * node.gaba * node.h_soil * ρ_H₂O(FT);
    m_out  = flow * Δt;
    m_all -= m_out;
    swc    = m_all / (node.gaba * node.h_soil * ρ_H₂O(FT));
    swc    = max(swc, FT(1e-6));

    soil_moisture_swc!(node, swc);

    return nothing
end
