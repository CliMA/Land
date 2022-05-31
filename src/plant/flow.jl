
abstract type AbstractFlowProfile{FT<:AbstractFloat} end



mutable struct NonSteadyStateFlow{FT} <: AbstractFlowProfile{FT}
    "Vector of buffer water flow `[mol m⁻²]`"
    f_buffer::Vector{FT}
    "Vector of xylem water flow `[mol m⁻²]`"
    f_element::Vector{FT}
    "Flow rate in `[mol s⁻¹]` or `[mol m⁻² s⁻¹]` (for leaf)"
    f_in::FT
    "Flow rate out `[mol s⁻¹]` or `[mol m⁻² s⁻¹]` (for leaf)"
    f_out::FT
end



NonSteadyStateFlow{FT}(N::Int, isleaf::Bool = true) where {FT<:AbstractFloat} = (
    if isleaf
        _zeros = zeros(FT,1);
    else
        _zeros = zeros(FT,N);
    end;

    return NonSteadyStateFlow{FT}(
                _zeros,     # f_buffer
                _zeros,     # f_element
                0,          # f_in
                0           # f_out
    )
);



mutable struct SteadyStateFlow{FT} <: AbstractFlowProfile{FT}
    "Flow rate `[mol s⁻¹]` or `[mol m⁻² s⁻¹]` (for leaf)"
    flow::FT
end
