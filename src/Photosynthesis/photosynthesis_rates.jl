abstract type AbstractPhotosynthesisModel end

abstract type AbstractC3Photosynthesis <: AbstractPhotosynthesisModel end
abstract type AbstractC4Photosynthesis <: AbstractPhotosynthesisModel end

struct  C3FvCBPhoto    <: AbstractC3Photosynthesis end
struct  C3FvCBPhotoGs  <: AbstractC3Photosynthesis end
struct  C3FvCBPhotoATP <: AbstractC3Photosynthesis end
struct  C4CollatzPhoto <: AbstractC4Photosynthesis end

"Compute ETR (Je) from Jmax and Jpot using quadratic"
function electron_transport_rate!(model::AbstractPhotosynthesisModel, flux, leaf)
    @unpack  Γstar, Cc, Vcmax, maxPSII, Jmax,θ_j, PSII_frac = leaf
    Je_pot = PSII_frac * maxPSII * flux.APAR;
    leaf.Je = lower_quadratic(θ_j, -(Je_pot + Jmax), Je_pot * Jmax)
end

############  model::C3FvCBPhoto #################################
"""
    rubisco_limited_rate(f::C3FvCBPhoto, leaf)
Solution when Rubisco activity is limiting (Classical Farquhar, von Caemmerer, Berry model) 
"""
function rubisco_limited_rate(model::AbstractC3Photosynthesis, flux, leaf)
    @unpack  Γstar, Cc, Kc, Ko, o₂, Vcmax = leaf
    Vcmax * max(Cc-Γstar, 0) / (Cc + Kc*(1 +o₂/Ko))
end

"""
    product_limited_rate(f::C3FvCBPhoto, leaf)
Solution when Triose Phosphate Export is limiting (CLM implementation)) 
"""
function product_limited_rate(model::AbstractC3Photosynthesis, leaf)
    leaf.Vcmax/2
end

"""
    light_limited_rate(f::C3FvCBPhoto, flux, leaf)
Solution when Electron Transport Rate is limiting (Classical Farquhar, von Caemmerer, Berry model)
C3: RuBP-limited photosynthesis (this is the NADPH requirement stochiometry)
Using curvature θ_j and Jmax
"""
function light_limited_rate(model::AbstractC3Photosynthesis, flux, leaf)
    @unpack  Γstar, Cc = leaf
    electron_transport_rate!(model, flux, leaf) 
    leaf.Je * max(Cc-Γstar, 0) / (4Cc + 8Γstar)
end

############  model::C3FvCBPhotoATP #################################
"""
    light_limited_rate(f::C3FvCBPhotoATP, leaf)
Solution when Electron Transport Rate is limiting (Classical Farquhar, von Caemmerer, Berry model)
C3: RuBP-limited photosynthesis (this is the ATP requirement stochiometry)
Using curvature θ_j and Jmax
"""
function light_limited_rate(model::C3FvCBPhotoATP, flux, leaf)
    @unpack  Γstar, Cc = leaf
    electron_transport_rate!(model, flux, leaf)
    leaf.Je * max(Cc-Γstar, 0) / (4.5*Cc + 10.5*Γstar)
end


############  model::C3FvCBPhotoGs #################################
"""
    rubisco_limited_rate(f::C3FvCBPhotoGs, leaf)
Solution when Rubisco activity is limiting and a given gs
Solves quadratic equation 
"""
function rubisco_limited_rate(model::C3FvCBPhotoGs, flux, leaf)
    @unpack  Γstar, Cc, Kc, Ko, o₂, Vcmax, gs, gm, Rd, Γstar = leaf
    @unpack  Ca,ra,g_m_s_to_mol_m2_s = flux
    d = Kc*(1+o₂/Ko)
    a = (ra/g_m_s_to_mol_m2_s + 1.6/gs + 1.0/gm) # = 1/gleaf
    b = -(Ca + d) - (Vcmax-Rd)*a
    c = Vcmax * (Ca - Γstar) - Rd * (Ca+d)
    lower_quadratic(a,b,c)
end

"""
    light_limited_rate(f::C3FvCBPhotoGs, flux, leaf)
Solution when Electron Transport Rate is limiting and a given gs
C3: RuBP-limited photosynthesis (this is the NADPH requirement stochiometry)
Using curvature θ_j and Jmax
"""
function light_limited_rate(model::C3FvCBPhotoGs, flux, leaf)
    @unpack  Γstar, Cc, gs, gm, Rd = leaf
    @unpack  Ca, ra, g_m_s_to_mol_m2_s = flux
    electron_transport_rate!(model, flux, leaf)
    a = 4*(ra/g_m_s_to_mol_m2_s + 1.6/gs + 1.0/gm) # = 1/gleaf
    b = -(4Ca + 8Γstar) - (leaf.Je - 4Rd)*a/4
    c = leaf.Je * (Ca-Γstar) - Rd * (4Ca + 8Γstar)
    lower_quadratic(a,b,c)
end

############  model::C4_Collatz_Photo #################################
"""
    rubisco_limited_rate(f::C4_Collatz_Photo, leaf)
Solution when C4 Rubisco activity is limiting (Collatz) 
"""
function rubisco_limited_rate(f::C4CollatzPhoto, flux, leaf)
    leaf.Vcmax
end

"""
    light_limited_rate(f::C3FvCBPhoto, leaf)
Solution when C4 Electron Transport Rate is limiting (Collatz model)
Linear with J
"""
function light_limited_rate(model::C4CollatzPhoto, flux, leaf)
    FT = typeof(APAR)
    # Efficiency of ETR for C4 plants ()
    α= FT(1/6)
    Je = α*(leaf.maxPSII * flux.APAR)/2;
end

"""
    product_limited_rate(f::C4_Collatz_Photo, leaf)
Solution for PEP carboxylase-limited rate of carboxylation
"""
function product_limited_rate(model::C4CollatzPhoto, leaf)
    @unpack  Kp,Vpmax, Cc = leaf
    Vpmax * CC / (Cc + Kp)
    # Still to be implemented, needs k_p * C_i/P
    Inf
end











