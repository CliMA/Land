abstract type AbstractFluorescenceModel end

Base.@kwdef struct FlexasTolBerryFluorescence{FT} <: AbstractFluorescenceModel 
    Kn1::FT = 5.01
    Kn2::FT = 1.93
    Kn3::FT = 10.0
end

"""
    leaf_fluorescence!(model,ps,x,leaf::leaf_params )

Compute Fluorescence yields, Kn and Kp.

# Arguments
- `model`: Model Type
- `ps::Number`: PSII yield.
- `x::Number`: Degree of light saturation: [0-1] .
- `leaf::leaf_params`: leaf_params structure.
"""
function leaf_fluorescence!(model::FlexasTolBerryFluorescence, ps::Number, x::Number, leaf::leaf_params)
    FT = eltype(ps)
    @unpack Kf,Kn,Kd = leaf
    @unpack Kn1,Kn2,Kn3 = model
    # Max PSII rate constant
    Kp_max = FT(4.0)

    x_alpha = exp(log(x)*Kn2);
    #println(x_alpha)
    leaf.Kn   = Kn1 * (1+Kn3)* x_alpha/(Kn3 + x_alpha);
    leaf.Kp   = max(0,-ps*(Kf+Kd+Kn)/(ps-1));
    leaf.Fo   = Kf/(Kf+Kp_max+Kd);
    leaf.Fo′  = Kf/(Kf+Kp_max+Kd+Kn);
    leaf.Fm   = Kf/(Kf   +Kd);
    leaf.Fm′  = Kf/(Kf   +Kd+Kn);
    leaf.ϕs   = leaf.Fm′*(1-ps);
    # leaf.eta  = leaf.ϕs/leaf.Fo; # don't need this anymore, better to use ϕs directly for SIF as Fo is not always fqe=0.01.
    leaf.qQ   = 1-(leaf.ϕs-leaf.Fo′)/(leaf.Fm-leaf.Fo′);
    leaf.qE   = 1-(leaf.Fm-leaf.Fo′)/(leaf.Fm′-leaf.Fo);

    leaf.NPQ  = Kn/(Kf+Kd);
end
