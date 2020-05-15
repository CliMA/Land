abstract type AbstractFluorescenceModel end

Base.@kwdef struct FlexasTolBerryFluorescence{FT} <: AbstractFluorescenceModel 
    Kn1::FT = 5.01
    Kn2::FT = 1.93
    Kn3::FT = 10.0
end

"""
    leaf_fluorescence!(model,leaf::leaf_params )

Compute Fluorescence yields, Kn and Kp.

# Arguments
- `model`: Model Type
- `leaf::leaf_params`: leaf_params structure.
"""
function leaf_fluorescence!(model::FlexasTolBerryFluorescence, leaf::leaf_params)
    #FT = eltype(ps)
    @unpack Kf,Kn,Kd = leaf
    @unpack Kn1,Kn2,Kn3 = model
    
    leaf.CO2_per_electron = (leaf.Cc-leaf.Γstar)/(leaf.Cc+2leaf.Γstar) * leaf.effcon;
    # Actual effective ETR:
    leaf.Ja = max(0,leaf.Ag / leaf.CO2_per_electron);
    leaf.Ja = min(leaf.Ja,leaf.Je_pot )
    # Effective photochemical yield:
    leaf.φ = leaf.maxPSII*leaf.Ja/leaf.Je_pot;

    #println(flux.Ja, " ", flux.Je_pot)
    leaf.φ = min(1/leaf.maxPSII,leaf.φ)

    x   = max(0,  1-leaf.φ/leaf.maxPSII);       # degree of light saturation: 'x' (van der Tol e.Ap. 2014)

    # Max PSII rate constant
    Kp_max = FT(4.0)

    x_alpha = exp(log(x)*Kn2);
    #println(x_alpha)
    leaf.Kn   = Kn1 * (1+Kn3)* x_alpha/(Kn3 + x_alpha);
    leaf.Kp   = max(0,-leaf.φ*(Kf+Kd+Kn)/(leaf.φ-1));

    leaf.Fo   = Kf/(Kf+Kp_max+Kd   );
    leaf.Fo′  = Kf/(Kf+Kp_max+Kd+Kn);
    leaf.Fm   = Kf/(Kf       +Kd   );
    leaf.Fm′  = Kf/(Kf       +Kd+Kn);
    leaf.ϕs   = leaf.Fm′*(1-leaf.φ);
    # leaf.eta  = leaf.ϕs/leaf.Fo; # don't need this anymore, better to use ϕs directly for SIF as Fo is not always fqe=0.01.
    leaf.qQ   = 1-(leaf.ϕs-leaf.Fo′)/(leaf.Fm-leaf.Fo′);
    leaf.qE   = 1-(leaf.Fm-leaf.Fo′)/(leaf.Fm′-leaf.Fo);

    leaf.NPQ  = Kn/(Kf+Kd);
end
