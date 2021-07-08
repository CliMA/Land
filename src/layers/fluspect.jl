###############################################################################
#
# Calculate leaf optical properties
# TODO need to be refactored for better documentation
#
###############################################################################
"""
    fluspect!(leaf::LeafBios{FT},
              wls::WaveLengths{FT};
              APAR_car::Bool = true
    ) where {FT<:AbstractFloat}

Computes leaf optical properties (reflectance and transittance) based on
    pigment concentrations. Also computes Fluorescence excitation matrices.
    Mostly based on PROSPECT-D for leaf reflectance/transmission and FluSpec
    for fluorescence.
- `leaf` [`LeafBios`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
- `APAR_car` If true, include Car absorption in APAR for photosynthesis
"""
function fluspect!(
            leaf::LeafBios{FT},
            wls::WaveLengths{FT};
            APAR_car::Bool = true
) where {FT<:AbstractFloat}
    # ***********************************************************************
    # Jacquemoud S., Baret F. (1990), PROSPECT: a model of leaf optical
    # properties spectra; Remote Sens. Environ.; 34:75-91.
    # Reference:
    # Féret, Gitelson, Noble & Jacquemoud [2017]. PROSPECT-D: Towards modeling
    # leaf optical properties through a complete lifecycle
    # Remote Sensing of Environment; 193:204215
    # DOI: http://doi.org/10.1016/j.rse.2017.03.004
    # The specific absorption coefficient corresponding to brown pigment is()
    # provided by Frederic Baret [EMMAH, INRA Avignon, baret@avignon.inra.fr]
    # & used with his autorization.
    # ***********************************************************************

    @unpack N, Cab, Car, Ant, Cs, Cw, Cm, ρ_SW, τ_SW, Cx, ndub = leaf;
    @unpack iWLE, iWLF, optis, WLE, WLF = wls;
    @unpack Kab, Kant, KBrown, Kcar, KcaV, KcaZ, Km, Kw, nr, phi = optis;

    #println(N, " ", Cab, " ", Car," ",  Ant, " ", Cs, " ", Cw, " ", Cm)
    Kcaro = (1 -Cx)* KcaV + Cx * KcaZ;

    Kall    = (Cab*Kab.+Car*Kcar.+Ant*Kant.+Cs*KBrown.+Cw*Kw.+Cm*Km)/N

    # Relative absorption by Chlorophyll and Carotenoids only (drives SIF and GPP eventually)
    if APAR_car
        leaf.kChlrel = (Cab*Kab+Car*Kcar)./(Kall*N.+eps(FT));
    else
        leaf.kChlrel = (Cab*Kab)./(Kall*N.+eps(FT));
    end
    leaf.kChlrel_old  = (Cab*Kab)./(Kall*N.+eps(FT));
    #println(typeof(Kall))
    # Adding eps() here to keep it stable and NOT set to 1 manually when Kall=0 (ForwardDiff won't work otherwise)
    tau = (1 .-Kall).*exp.(-Kall) .+ Kall.^2 .*real.(expint.(Kall.+eps(FT)))

    # ***********************************************************************
    # reflectance & transmittance of one layer
    # ***********************************************************************
    # Allen W.A., Gausman H.W., Richardson A.J., Thomas J.R. (1969)
    # Interaction of isotropic ligth with a compact plant leaf; J. Opt.
    # Soc. Am., 59[10]:1376-1379.
    # ***********************************************************************
    # reflectivity & transmissivity at the interface
    #-------------------------------------------------
    # From Prospect-D, uses 40 here instead of 59 from CVT)

    #talf    = calctav.(59.,nr)
    talf    = calctav.(FT(40),nr)
    ralf    = FT(1) .-talf

    t12     = calctav.(FT(90), nr)
    r12     = FT(1) .-t12
    t21     = t12./(nr.^2)
    r21     = FT(1) .-t21
    #println(typeof(r21))
    # top surface side
    denom   = FT(1) .-r21.*r21.*tau.^2
    Ta      = talf.*tau.*t21./denom
    Ra      = ralf.+r21.*tau.*Ta
    #println(typeof(t12), typeof(t21), typeof(denom), typeof(tau))
    # bottom surface side
    t       = t12.*tau.*t21./denom
    r       = r12+r21.*tau.*t

    # ***********************************************************************
    # reflectance & transmittance of N layers
    # Stokes equations to compute properties of next N-1 layers [N real]
    # Normal case()
    # ***********************************************************************
    # Stokes G.G. (1862), On the intensity of the light reflected from
    # | transmitted through a pile of plates; Proc. Roy. Soc. Lond.
    # 11:545-556.
    # ***********************************************************************
    D       = sqrt.((FT(1) .+r.+t).*(FT(1) .+r.-t).*(FT(1) .-r.+t).*(FT(1) .-r.-t))
    #println(typeof(D), typeof(r), typeof(t))
    rq      = r.^2
    tq      = t.^2
    a       = (FT(1) .+rq.-tq.+D)./(2r)
    b       = (FT(1) .-rq.+tq.+D)./(2t)

    bNm1    = b.^(N-1);                  #
    bN2     = bNm1.^2
    a2      = a.^2
    denom   = a2.*bN2.-1
    Rsub    = a.*(bN2.-1)./denom
    Tsub    = bNm1.*(a2.-1)./denom

    # Case of zero absorption
    j       = findall(r.+t .>= 1)
    Tsub[j] = t[j]./(t[j]+(1 .-t[j])*(leaf.N-1))
    Rsub[j] = 1 .-Tsub[j]

    # Reflectance & transmittance of the leaf: combine top layer with next N-1 layers
    denom   = 1 .-Rsub.*r
    leaf.τ_SW    = Ta.*Tsub./denom
    leaf.ρ_SW    = Ra.+Ta.*Rsub.*t./denom
    leaf.α_SW    = 1 .- leaf.τ_SW .- leaf.ρ_SW;
    τ_SW    = leaf.τ_SW
    ρ_SW    = leaf.ρ_SW
    #RT     = [refl tran]
    #
    #
    #
    #
    # a place holder
    #if leaf.fqe ==0.0
    #    # do something?
    #end
    #
    #
    #
    #
    # FROM SCOPE notes:
    # From here a new path is taken: The doubling method used to calculate
    # fluoresence is now only applied to the part of the leaf where absorption
    # takes place, that is, the part exclusive of the leaf-air interfaces. The
    # reflectance (rho) and transmittance (tau) of this part of the leaf are
    # now determined by "subtracting" the interfaces
    # CF Note: All of the below takes about 10 times more time than the RT above. Need to rething speed and accuracy. (10nm is bringing it down a lot!)

    Rb  = (ρ_SW-ralf)./(talf.*t21+(ρ_SW-ralf).*r21);  # Remove the top interface
    tt1 = (talf.*t21);
    tt2 = τ_SW.*(1 .-Rb.*r21)
    Z   = tt2./tt1;            # Derive Z from the transmittance

    tt1 = (Rb-r21.*Z.^2); tt2 = (1 .-(r21.*Z).^2)
    rho = tt1./tt2;      # Reflectance and transmittance
    tt1 = (1 .-Rb.*r21);
    tau = tt1./tt2.*Z;    # of the leaf mesophyll layer
    t   =   tau;
    r   =   max.(rho,0);                       # Avoid negative r

    # Derive Kubelka-Munk s and k
    I_rt     =   findall((r.+t).<1);
    D[I_rt]  =   sqrt.((1 .+ r[I_rt] .+ t[I_rt]) .* (1 .+ r[I_rt] .- t[I_rt]) .* (1 .- r[I_rt] .+ t[I_rt]) .*  (1 .- r[I_rt] .- t[I_rt]));
    a[I_rt]  =   (1 .+ r[I_rt].^2 .- t[I_rt].^2 .+ D[I_rt]) ./ (2r[I_rt]);
    b[I_rt]  =   (1 .- r[I_rt].^2 + t[I_rt].^2 .+ D[I_rt]) ./ (2t[I_rt]);
    a[(r.+t).>=1] .=   FT(1.0);
    b[(r.+t).>=1] .=   FT(1.0);

    s        =   r./t;
    I_a      =   findall((a.>1).&(a.!=Inf));
    s[I_a]   =   2 .*a[I_a] ./ (a[I_a].^2 .- 1) .* log.(b[I_a]);

    k        =   log.(b);
    k[I_a]   =   (a[I_a].-1) ./ (a[I_a].+1) .* log.(b[I_a]);
    kChl     =   leaf.kChlrel .* k;

    # indices of WLE and WLF within wlp

    epsi         = FT(2)^(-ndub);

    # initialisations
    te     = 1 .-(k[iWLE].+s[iWLE]) * epsi;
    tf     = 1 .-(k[iWLF].+s[iWLF]) * epsi;
    re     = s[iWLE] * epsi;
    rf     = s[iWLF] * epsi;

    sigmoid     = 1 ./(1 .+exp.(-WLF/10).*exp.(WLE'/10));  # matrix computed as an outproduct
    #println(size(sigmoid)," ", size(phi), " ", size(kChl)," ", size(iWLE), " ", size(iWLF), " ", size(kChl[iWLE]))

    Mf = leaf.fqe .* ((FT(0.5)*phi[iWLF]).*epsi) .* kChl[iWLE]'.*sigmoid
    Mb = leaf.fqe .* ((FT(0.5)*phi[iWLF]).*epsi) .* kChl[iWLE]'.*sigmoid

    Ih          = ones(FT,1,length(te));     # row of ones
    Iv          = ones(FT,length(tf),1);     # column of ones

    # Doubling Adding Routine
    for i = 1:ndub
        xe = te./(1 .-re.*re);  ten = te.*xe;  ren = re.*(1 .+ten);
        xf = tf./(1 .-rf.*rf);  tfn = tf.*xf;  rfn = rf.*(1 .+tfn);

        A11  = xf*Ih + Iv*xe';
        A12 = (xf*xe').*(rf*Ih .+ Iv*re');
        A21  = 1 .+(xf*xe').*(1 .+rf*re');
        A22 = (xf.*rf)*Ih+Iv*adjoint(xe.*re);
        #println(typeof(opti.phi), typeof(kChl), typeof(Mf))
        Mfn   = Mf  .* A11 .+ Mb  .* A12;
        Mbn   = Mb  .* A21 .+ Mf  .* A22;

        te   = ten;  re  = ren;   tf   = tfn;   rf   = rfn;
        Mf  = Mfn; Mb = Mbn;
    end
    # Here we add the leaf-air interfaces again for obtaining the final
    # leaf level fluorescences.

    # This heree reduced red SIF quite a bit in backscatter, not sure why.
    g = Mb; f = Mf;

    Rb = rho .+ tau.^2 .*r21./(1 .-rho.*r21);

    Xe = Iv * adjoint(talf[iWLE]./(1 .-r21[iWLE].*Rb[iWLE]));
    Xf = t21[iWLF]./(1 .-r21[iWLF].*Rb[iWLF]) * Ih;
    Ye = Iv * adjoint(tau[iWLE].*r21[iWLE]./(1 .-rho[iWLE].*r21[iWLE]));
    Yf = tau[iWLF].*r21[iWLF]./(1 .-rho[iWLF].*r21[iWLF]) * Ih;

    A = Xe .* (1 .+ Ye.*Yf) .* Xf;
    B = Xe .* (Ye .+ Yf) .* Xf;

    gn = A .* g + B .* f;
    fn = A .* f + B .* g;

    leaf.Mb  = gn;
    leaf.Mf  = fn;

    return nothing
end
