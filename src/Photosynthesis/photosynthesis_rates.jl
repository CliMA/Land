abstract type AbstractPhotosynthesisModel end

abstract type  C3_FvCB_Photo <: AbstractPhotosynthesisModel end
abstract type  C3_FvCB_Photo_gs <: C3_FvCB_Photo end
abstract type  C3_FvCB_PhotoATP <: C3_FvCB_Photo end
abstract type  C4_Collatz_Photo <: AbstractPhotosynthesisModel end

"Compute ETR (Je) from Jmax and Jpot using quadratic"
function electron_transport_rate!(model::AbstractPhotosynthesisModel, flux, leaf)
    @unpack  Γstar, Cc, Vcmax, maxPSII, Jmax,θ_j, PSII_frac = leaf
    Je_pot = leaf.PSII_frac * (leaf.maxPSII * flux.APAR);
    leaf.Je = lower_quadratic(θ_j, -(Je_pot + Jmax), Je_pot * Jmax)
end

############  model::C3_FvCB_Photo #################################
"""
    rubisco_limited_rate(f::C3_FvCB_Photo, leaf)
Solution when Rubisco activity is limiting (Classical Farquhar, von Caemmerer, Berry model) 
"""
function rubisco_limited_rate(model::C3_FvCB_Photo, flux, leaf)
    @unpack  Γstar, Cc, Kc, Ko, o₂, Vcmax = leaf
    Vcmax * max(Cc-Γstar, 0) / (Cc + Kc*(1 +o₂/Ko))
end

"""
    product_limited_rate(f::C3_FvCB_Photo, leaf)
Solution when Triose Phosphate Export is limiting (CLM implementation)) 
"""
function product_limited_rate(model::C3_FvCB_Photo, leaf)
    leaf.Vcmax/2
end

"""
    light_limited_rate(f::C3_FvCB_Photo, flux, leaf)
Solution when Electron Transport Rate is limiting (Classical Farquhar, von Caemmerer, Berry model)
C3: RuBP-limited photosynthesis (this is the NADPH requirement stochiometry)
Using curvature θ_j and Jmax
"""
function light_limited_rate(model::C3_FvCB_Photo, flux, leaf)
    @unpack  Γstar, Cc = leaf
    electron_transport_rate(model, flux, leaf) * max(Cc-Γstar, 0) / (4Cc + 8Γstar)
end

############  model::C3_FvCB_PhotoATP #################################
"""
    light_limited_rate(f::C3_FvCB_PhotoATP, leaf)
Solution when Electron Transport Rate is limiting (Classical Farquhar, von Caemmerer, Berry model)
C3: RuBP-limited photosynthesis (this is the ATP requirement stochiometry)
Using curvature θ_j and Jmax
"""
function light_limited_rate(model::C3_FvCB_PhotoATP, flux, leaf)
    @unpack  Γstar, Cc = leaf
    electron_transport_rate!(model, flux, leaf)
    leaf.Je * max(Cc-Γstar, 0) / (4.5*Cc + 10.5*Γstar)
end


############  model::C3_FvCB_Photo_gs #################################
"""
    rubisco_limited_rate(f::C3_FvCB_Photo_gs, leaf)
Solution when Rubisco activity is limiting and a given gs
Solves quadratic equation 
"""
function rubisco_limited_rate(model::C3_FvCB_Photo_gs, flux, leaf)
    @unpack  Γstar, Cc, Kc, Ko, o₂, Vcmax, gs, gm, Rd, Γstar = leaf
    @unpack  Ca,ra,g_m_s_to_mol_m2_s = flux
    d = Kc*(1+o₂/Ko)
    a = (ra/g_m_s_to_mol_m2_s + 1.6/gs + 1.0/gm) # = 1/gleaf
    b = -(Ca + d) - (Vcmax-Rd)*a
    c = Vcmax * (Ca - Γstar) - Rd * (Ca+d)
    lower_quadratic(a,b,c)
end

"""
    light_limited_rate(f::C3_FvCB_Photo_gs, flux, leaf)
Solution when Electron Transport Rate is limiting and a given gs
C3: RuBP-limited photosynthesis (this is the NADPH requirement stochiometry)
Using curvature θ_j and Jmax
"""
function light_limited_rate(model::C3_FvCB_Photo_gs, flux, leaf)
    @unpack  Γstar, Cc = leaf
    @unpack  Ca,ra,g_m_s_to_mol_m2_s = flux
    electron_transport_rate!(model, flux, leaf) * max(Cc-Γstar, 0) / (4Cc + 8Γstar)
    a = (ra/g_m_s_to_mol_m2_s + 1.6/gs + 1.0/gm) # = 1/gleaf
    b = -(4Ca + 8Γstar) - (Je - 4Rd)*a
    c = leaf.Je * (Ca-Γstar) - Rd * (4Ca + 8Γstar)
    lower_quadratic(a,b,c)
end

############  model::C4_Collatz_Photo #################################
"""
    rubisco_limited_rate(f::C4_Collatz_Photo, leaf)
Solution when C4 Rubisco activity is limiting (Collatz) 
"""
function rubisco_limited_rate(f::C4_Collatz_Photo, flux, leaf)
    leaf.Vcmax
end

"""
    light_limited_rate(f::C3_FvCB_Photo, leaf)
Solution when C4 Electron Transport Rate is limiting (Collatz model)
Linear with J
"""
function light_limited_rate(model::C4_Collatz_Photo, flux, leaf)
    FT = typeof(APAR)
    # Efficiency of ETR for C4 plants ()
    α= FT(1/6)
    Je = α*(leaf.maxPSII * flux.APAR)/2;
end

"""
    product_limited_rate(f::C4_Collatz_Photo, leaf)
Solution for PEP carboxylase-limited rate of carboxylation
"""
function product_limited_rate(model::C4_Collatz_Photo, leaf)
    @unpack  Kp,Vpmax, Cc = leaf
    Vpmax * CC / (Cc + Kp)
    # Still to be implemented, needs k_p * C_i/P
    Inf
end











