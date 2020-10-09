###############################################################################
#
# run ODE on vertical layered trace gasses
#
###############################################################################
"""
    vertical_layers!(vls::VerticalLayers{FT}, t::FT) where {FT<:AbstractFloat}

Update vertical trace gas information with time, given
- `vls` [`VerticalLayers`](@ref) type struct
- `t` Time since last update
"""
function vertical_layers!(
            vls::VerticalLayers{FT},
            t::FT
) where {FT<:AbstractFloat}
    # unpack parameters
    @unpack d_CO₂, d_layer, p_CO₂, x_layer, Δ_CO₂ = vls;
    bc = DirichletBC(FT(50), FT(30));
    u0 = vls.p_CO₂;

    # define step function, move out later
    # very memory extensive for now
    step(u,p,t) = Δ_CO₂ * bc * u .* d_CO₂ ./ x_layer ./ d_layer;
    prob = ODEProblem(step, u0, (FT(0), t));
    alg  = KenCarp4();
    sol  = solve(prob, alg);

    @show sol[end];

    return nothing
end
