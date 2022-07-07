
function ∂g∂t end;



∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}, β_param::BetaParameterG1; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaf, air; β = β);

    return (_gsw - leaf.g_H₂O_s) / sm.Τ
);



∂g∂t(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}, leaf::Leaf{FT}, air::AirLayer{FT}, β_param::BetaParameterVcmax; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _gsw = empirical_equation(sm, leaf, air; β = FT(1));

    return (_gsw - leaf.g_H₂O_s) / sm.Τ
);





function stomatal_conductance! end



stomatal_conductance!(leaf::Leaf{FT}, air::AirLayer{FT}, Δt::FT = 60; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack SM = leaf;

    leaf.g_H₂O_s += ∂g∂t(SM, leaf, air, SM.Β; β = β) * Δt;

    return nothing
);
