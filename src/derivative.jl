#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jul-08: migrate function from StomataModels.jl
#     2022-Jul-08: add method for LeafHydraulics
#     2022-Jul-08: add methods for Leaf, Leaves1D, and Leaves2D
#
#######################################################################################################################################################################################################
"""

    ∂E∂P(lf::Union{Leaf{FT}, Leaves2D{FT}}, flow::FT; pm::Int = 1) where {FT<:AbstractFloat}
    ∂E∂P(lf::Leaves1D{FT}, flow::FT, ind::Int; pm::Int = 1) where {FT<:AbstractFloat}

Return the marginal hydraulic conductance, given
- `lf` `Leaf`, `Leaves1D`, or `Leaves2D` type struct
- `flow` Flow rate through the leaf xylem `[mol m⁻² s⁻¹]`
- `pm` Scaling factor for incremental flow rate, default is 1
- `ind` Which leaf in `Leaves1D` (1 for sunlit and 2 for shaded)

"""
function ∂E∂P end

∂E∂P(lf::Union{Leaf{FT}, Leaves2D{FT}}, flow::FT; pm::Int = 1) where {FT<:AbstractFloat} = ∂E∂P(lf.HS, flow, lf.t; pm = pm);

∂E∂P(lf::Leaves1D{FT}, flow::FT, ind::Int; pm::Int = 1) where {FT<:AbstractFloat} = ∂E∂P(lf.HS, flow, lf.t[ind]; pm = pm);

∂E∂P(hs::LeafHydraulics{FT}, flow::FT, t::FT; pm::Int = 1) where {FT<:AbstractFloat} = (
    @assert pm != 0 "pm must not be 0";

    # calculate the flow and pressure at the given flow rate
    _e1 = flow;
    _p1 = xylem_end_pressure(hs, _e1, t);

    # compute the P and E when e de- or in-creases by 1e-7 mol m⁻² s⁻¹
    _e2 = _e1 + pm * FT(1e-7);
    _p2 = xylem_end_pressure(hs, _e2, t);

    return -1 * (_e2 - _e1) / (_p2 - _p1)
);
