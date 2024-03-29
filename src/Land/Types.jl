#=
###############################################################################
#
# Layer of soil and air to run ODE
#
###############################################################################
"""
    mutable struct VerticalLayers{FT}

Struct that store trace gas information along vertical layers.

"""
Base.@kwdef mutable struct VerticalLayers{FT}
    "Number of layers"
    n_layer::Int = 20
    "Thickness per layer"
    d_layer::Array{FT,1} = FT[1 for i in 1:n_layer]
    "Vertical distances among layers (0.5 for boundaries)"
    x_layer::Array{FT,1} = FT[0.5; [1 for i in 1:n_layer-2]; 0.5]

    # Temperature
    "Temperature at each layer"
    T::Array{FT,1} = ones(FT,n_layer) .* FT(298.15)

    # CO₂ partial pressure
    "Diffusion coefficient array"
    d_CO₂::Vector{FT} = diffusive_coefficient.(T, [TraceGasCO₂{FT}()], [TraceGasAir{FT}()])
    "Vertical CO₂ partial pressure"
    p_CO₂::Vector{FT} = ones(FT,n_layer) .* 41
    "derivative operator for CO₂ partial pressure"
    Δ_CO₂::DerivativeOperator = CenteredDifference(2, 2, FT(1), n_layer)
end
=#
