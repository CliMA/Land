#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jul-08: migrate function from StomataModels.jl
#     2022-Jul-08: add method for LeafHydraulics
#     2022-Jul-08: add methods for Leaf, Leaves1D, and Leaves2D
#     2022-Jul-08: add option δe to be more general
#
#######################################################################################################################################################################################################
"""

    ∂E∂P(lf::Union{Leaf{FT}, Leaves2D{FT}}, flow::FT; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
    ∂E∂P(lf::Leaves1D{FT}, flow::FT, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}

Return the marginal hydraulic conductance, given
- `lf` `Leaf`, `Leaves1D`, or `Leaves2D` type struct
- `flow` Flow rate through the leaf xylem `[mol m⁻² s⁻¹]`
- `δe` Incremental flow rate, default is 1e-7
- `ind` Which leaf in `Leaves1D` (1 for sunlit and 2 for shaded)

"""
function ∂E∂P end

∂E∂P(lf::Union{Leaf{FT}, Leaves2D{FT}}, flow::FT; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = ∂E∂P(lf.HS, flow, lf.t; δe = δe);

∂E∂P(lf::Leaves1D{FT}, flow::FT, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = ∂E∂P(lf.HS, flow, lf.t[ind]; δe = δe);

∂E∂P(hs::LeafHydraulics{FT}, flow::FT, t::FT; δe::FT = FT(1e-7)) where {FT<:AbstractFloat} = (
    @assert δe != 0 "pm must not be 0";

    _p1 = xylem_end_pressure(hs, flow, t);
    _p2 = xylem_end_pressure(hs, flow + δe, t);

    return -δe / (_p2 - _p1)
);
