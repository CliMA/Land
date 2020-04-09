module CanopyRTMod
#using GSL
using Polynomials
using Statistics
# Matlab reading
using MAT
# Numerical integration package (Simpson rule)
using QuadGK
using MathToolsMod

export optipar, angles, leafbio, angle, leaf, optis

# Floating point precision can be changed here
FT = Float32

include("PhotoStructs.jl")

# Leaf inclination distribution
const litab   = FT[5.,15.,25.,35.,45.,55.,65.,75.,81.,83.,85.,87.,89.];

const minwlPAR  = FT(400.); # PAR range
const maxwlPAR  = FT(700.);

const minwle  = FT(400.); # SIF excitation range
const maxwle  = FT(750.);

const minwlf  = FT(650.); # SIF emission range
const maxwlf  = FT(850.);

const fac     = FT(0.001) #(unit conversion factor mW >> W)

# Doubling Adding layers
const ndub = 15
const nLayers = 20
# Load some standard settings here (can be done differently later but not bad to be fixed here)
swl = collect(FT(400):10:2401)

# Load optical properties (these can be fixed globally!):
optis = optipar{FT}(loadOpti(swl)...);

# Set wavelength here (will stay constant!!)
const wl = optis.lambda
const nwl = length(wl)
const Iwle   = findall((wl.>=minwle) .& (wl.<=maxwle));
const Iwlf   = findall((wl.>=minwlf) .& (wl.<=maxwlf));
const iPAR   = findall((wl.>=minwlPAR) .& (wl.<=maxwlPAR))
const wle    = wl[Iwle];    # excitation wavelengths,  column
const wlf    = wl[Iwlf];    # fluorescence wavelengths,  column

# Later on, we can implement an array of these structs for global parallel calculations!
# Load sun (will be coming from CLIMA in future, just file for now)
sunRad = incomingRadiation{FT}(loadSun(swl)...)
leaf   = leafbio{FT}(ρ_SW=similar(wl),τ_SW=similar(wl), kChlrel=similar(wl), Mb=zeros(FT,length(wle), length(wlf)), Mf=zeros(FT,length(wle), length(wlf)) )
angle  = struct_angles{FT}()
canopy = struct_canopy{FT}()
# Canopy radiation:
canRad = getRadStruct(length(wl), canopy.nlayers, length(canopy.lazitab), length(litab),FT)
# Canopy Optical Properties:
canOpt = getCanOptStruct(length(wl), canopy.nlayers, length(canopy.lazitab), length(litab),FT)
const nlayers = canopy.nlayers

function fluspect!(leaf::leafbio, opti::optipar)
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

    Kcaro = (1 -leaf.Cx).* opti.KcaV + leaf.Cx .* opti.KcaZ;

    Kall    = (leaf.Cab*opti.Kab.+leaf.Car*opti.Kcar.+leaf.Ant*opti.Kant.+leaf.Cs*opti.KBrown.+leaf.Cw*opti.Kw.+leaf.Cm*opti.Km)./leaf.N
    # Relative absorption by Chlorophyll only (drives SIF and GPP eventually)
    leaf.kChlrel  = leaf.Cab*opti.Kab./(Kall.*leaf.N.+eps(FT));
    #println(typeof(Kall))
    # Adding eps() here to keep it stable and NOT set to 1 manually when Kall=0 (ForwardDiff won't work otherwise)
    tau = (1 .-Kall).*exp.(-Kall) .+ Kall.^2 .*real.(expint.(Kall.+eps(FT)))
    #println(typeof(expint.(Kall)))
    #println(typeof(real.((1 .-Kall).*exp.(-Kall)) .+ Kall.^2), typeof(expint.(Kall)))
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
    talf    = calctav.(40,opti.nr)
    ralf    = 1 .-talf

    t12     = calctav.(90,opti.nr)
    r12     = 1 .-t12
    t21     = t12./(opti.nr.^2)
    r21     = 1 .-t21

    # top surface side
    denom   = 1 .-r21.*r21.*tau.^2
    Ta      = talf.*tau.*t21./denom
    Ra      = ralf.+r21.*tau.*Ta
    #println(typeof(t12), typeof(t21), typeof(denom), typeof(tau))
    # bottom surface side
    t       = t12.*tau.*t21./denom
    r       = r12.+r21.*tau.*t

    # ***********************************************************************
    # reflectance & transmittance of N layers
    # Stokes equations to compute properties of next N-1 layers [N real]
    # Normal case()
    # ***********************************************************************
    # Stokes G.G. (1862), On the intensity of the light reflected from
    # | transmitted through a pile of plates; Proc. Roy. Soc. Lond.
    # 11:545-556.
    # ***********************************************************************
    D       = sqrt.((1 .+r.+t).*(1 .+r.-t).*(1 .-r.+t).*(1 .-r.-t))
    #println(typeof(D), typeof(r), typeof(t))
    rq      = r.^2
    tq      = t.^2
    a       = (1 .+rq.-tq.+D)./(2r)
    b       = (1 .-rq.+tq.+D)./(2t)

    bNm1    = b.^(leaf.N.-1);                  #
    bN2     = bNm1.^2
    a2      = a.^2
    denom   = a2.*bN2.-1
    Rsub    = a.*(bN2.-1)./denom
    Tsub    = bNm1.*(a2.-1)./denom

    # Case of zero absorption
    j       = findall(r.+t .>= 1)
    Tsub[j] = t[j]./(t[j]+(1 .-t[j])*(leaf.N-1))
    Rsub[j]	= 1 .-Tsub[j]

    # Reflectance & transmittance of the leaf: combine top layer with next N-1 layers
    denom   = 1 .-Rsub.*r
    leaf.τ_SW    = Ta.*Tsub./denom
    leaf.ρ_SW    = Ra.+Ta.*Rsub.*t./denom

    #RT     = [refl tran]
    if leaf.fqe ==0.0
        return
    end
    # FROM SCOPE notes:
    # From here a new path is taken: The doubling method used to calculate
    # fluoresence is now only applied to the part of the leaf where absorption
    # takes place, that is, the part exclusive of the leaf-air interfaces. The
    # reflectance (rho) and transmittance (tau) of this part of the leaf are
    # now determined by "subtracting" the interfaces
    # CF Note: All of the below takes about 10 times more time than the RT above. Need to rething speed and accuracy. (10nm is bringing it down a lot!)

    Rb  = (leaf.ρ_SW.-ralf)./(talf.*t21+(leaf.ρ_SW.-ralf).*r21);  # Remove the top interface
    Z   = leaf.τ_SW.*(1 .-Rb.*r21)./(talf.*t21);            # Derive Z from the transmittance

    rho = (Rb.-r21.*Z.^2)./(1 .-(r21.*Z).^2);      # Reflectance and transmittance
    tau = (1 .-Rb.*r21)./(1 .-(r21.*Z).^2).*Z;    # of the leaf mesophyll layer
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

    # indices of wle and wlf within wlp

    epsi         = FT(2)^(-ndub);

    # initialisations
    te     = 1 .-(k[Iwle].+s[Iwle]) .* epsi;
    tf     = 1 .-(k[Iwlf].+s[Iwlf]) .* epsi;
    re     = s[Iwle] .* epsi;
    rf     = s[Iwlf] .* epsi;

    sigmoid     = 1 ./(1 .+exp.(-wlf./10).*exp.(wle'./10));  # matrix computed as an outproduct
    #println(size(sigmoid)," ", size(phi), " ", size(kChl)," ", size(Iwle), " ", size(Iwlf), " ", size(kChl[Iwle]))

    Mf = Mb = leaf.fqe .* ((FT(0.5)*opti.phi[Iwlf]).*epsi) .* kChl[Iwle]'.*sigmoid

    Ih          = ones(FT,1,length(te));     # row of ones
    Iv          = ones(FT,length(tf),1);     # column of ones
    #A11 = A12 = A21 = A22 = zeros(FT,length(tf),length(te));
    #println(length(tf),length(te))
    # Doubling routine
    for i = 1:ndub
        xe = te./(1 .-re.*re);  ten = te.*xe;  ren = re.*(1 .+ten);
        xf = tf./(1 .-rf.*rf);  tfn = tf.*xf;  rfn = rf.*(1 .+tfn);

        A11  = xf*Ih + Iv*xe';
        A12 = (xf*xe').*(rf*Ih .+ Iv*re');
        A21  = 1 .+(xf*xe').*(1 .+rf*re');
        A22 = (xf.*rf)*Ih+Iv*(xe.*re)';
        #println(typeof(opti.phi), typeof(kChl), typeof(Mf))
        Mfn   = Mf  .* A11 .+ Mb  .* A12;
        Mbn   = Mb  .* A21 .+ Mf  .* A22;

        te   = ten;  re  = ren;   tf   = tfn;   rf   = rfn;
        Mf  = Mfn; Mb = Mbn;
    end

    # Here we add the leaf-air interfaces again for obtaining the final
    # leaf level fluorescences.
    g = Mb; f = Mf;

    Rb = rho .+ tau.^2 .*r21./(1 .-rho.*r21);

    Xe = Iv * (talf[Iwle]./(1 .-r21[Iwle].*Rb[Iwle]))';
    Xf = t21[Iwlf]./(1 .-r21[Iwlf].*Rb[Iwlf]) * Ih;
    Ye = Iv * (tau[Iwle].*r21[Iwle]./(1 .-rho[Iwle].*r21[Iwle]))';
    Yf = tau[Iwlf].*r21[Iwlf]./(1 .-rho[Iwlf].*r21[Iwlf]) * Ih;

    A = Xe .* (1 .+ Ye.*Yf) .* Xf;
    B = Xe .* (Ye .+ Yf) .* Xf;

    gn = A .* g + B .* f;
    fn = A .* f + B .* g;

    leaf.Mb  = gn;
    leaf.Mf  = fn;
    return
end

function RTM_SW!(can::struct_canopy, cO::struct_canopyOptProps, cR::struct_canopyRadiation, sun::incomingRadiation)
    nl = can.nlayers
    # Leaf Area Index
    LAI = can.LAI
    iLAI    = LAI/nl;
    # Define soil as polynomial (depends on state vector size), still TBD in the structure mode now.
    pSoil =  Polynomials.Poly(FT(0.2))
    rsoil = Polynomials.polyval(pSoil,wl.-mean(wl));
    soil_emissivity = 1 .-rsoil # emissivity of soil (goes into structure later)


    # 2.2. convert scattering and extinction coefficients into thin layer reflectances and transmittances
    # Eq. 17 in mSCOPE paper (changed here to compute real transmission)
    τ_ss = FT(exp(-cO.ks*iLAI));
    # Changed from mSCOPE paper here as well to exp. form but might be correct here.
    τ_dd = exp.(-cO.a.*iLAI);
    τ_sd = cO.sf.*iLAI;
    ρ_sd = cO.sb.*iLAI;
    ρ_dd = cO.sigb.*iLAI;

    # 2.3. reflectance calculation
    # Eq. 18 in mSCOPE paper

    Xss  = τ_ss;

    #Es_    = zeros(FT,nl+1,nwl)

    # Soil reflectance boundary condition (same for diffuse and direct)
    cO.R_sd[:,end] = rsoil
    cO.R_dd[:,end] = rsoil

    # Eq. 18 in mSCOPE paper
    @inbounds for j=nl:-1:1
        dnorm       = 1 .-ρ_dd[:,j].*cO.R_dd[:,j+1];
        cO.Xsd[:,j]    = (τ_sd[:,j]+τ_ss.*cO.R_sd[:,j+1].*ρ_dd[:,j])./dnorm;
        cO.Xdd[:,j]    = τ_dd[:,j]./dnorm;
        cO.R_sd[:,j]   = ρ_sd[:,j]+τ_dd[:,j].*(cO.R_sd[:,j+1].*Xss+cO.R_dd[:,j+1].*cO.Xsd[:,j]);
        cO.R_dd[:,j]   = ρ_dd[:,j]+τ_dd[:,j].*cO.R_dd[:,j+1].*cO.Xdd[:,j];
    end

    # 3.2 flux profile calculation
    # Eq. 19 in mSCOPE paper
    # Boundary condition at top: Incoming solar radiation
    cO.Es_[:,1]       = sun.E_direct;
    cR.E_down[:,1] = sun.E_diffuse;

    @inbounds for j=1:nl # from top to bottom
        cR.netSW_direct[:,j] = cO.Es_[:,j].*(1 .-(τ_ss.+τ_sd[:,j]+ρ_sd[:,j]));
        cO.Es_[:,j+1]    =   Xss.*cO.Es_[:,j];
        cR.E_down[:,j+1]  =   cO.Xsd[:,j].*cO.Es_[:,j]+cO.Xdd[:,j].*cR.E_down[:,j];
        cR.E_up[:,j]    =   cO.R_sd[:,j].*cO.Es_[:,j]+cO.R_dd[:,j].*cR.E_down[:,j];
    end
    # Boundary condition at the bottom, soil reflectance (Lambertian here)
    cR.E_up[:,end] = cO.R_sd[:,end].*cO.Es_[:,end]+cO.R_dd[:,end].*cR.E_down[:,end];

    # Hemispheric total outgoing (needs to go into stucture later):
    cR.Eout[:] = cR.E_up[:,1]

    # compute net diffuse radiation per layer:
    @inbounds for j=1:nl
        E_ = cR.E_down[:,j] +  cR.E_up[:,j+1];
        cR.netSW_diffuse[:,j]= E_.*(1 .-(τ_dd[:,j]+ρ_dd[:,j]))
    end
end


function computeCanopyMatrices!(leaf::Array{leafbio{FT}, 1},cO::struct_canopyOptProps)
    # 2. Calculation of reflectance
    # 2.1  reflectance, transmittance factors in a thin layer the following are vectors with length [nl,nwl]
    @inbounds for i=1:nlayers
        τ_SW = leaf[i].τ_SW
        ρ_SW = leaf[i].ρ_SW
        @inbounds for j=1:size(cO.sigb, 1)
          cO.sigb[j,i] = cO.ddb[j]  * ρ_SW[j]   + cO.ddf[j] * τ_SW[j]; # [nl,nwl]  diffuse     backscatter scattering coefficient for diffuse  incidence
          cO.sigf[j,i] = cO.ddf[j]  * ρ_SW[j]   + cO.ddb[j] * τ_SW[j]; # [nl,nwl]  diffuse     forward     scattering coefficient for diffuse  incidence
          cO.sb[j,i]   = cO.sdb[j]  * ρ_SW[j]   + cO.sdf[j] *  τ_SW[j]; # [nl,nwl]  diffuse     backscatter scattering coefficient for specular incidence
          cO.sf[j,i]   = cO.sdf[j]  * ρ_SW[j]   + cO.sdb[j] *  τ_SW[j]; # [nl,nwl]  diffuse     forward     scattering coefficient for specular incidence
          cO.vb[j,i]   = cO.dob[j]  * ρ_SW[j]   + cO.dof[j] *  τ_SW[j]; # [nl,nwl]  directional backscatter scattering coefficient for diffuse  incidence
          cO.vf[j,i]   = cO.dof[j]  * ρ_SW[j]   + cO.dob[j] *  τ_SW[j]; # [nl,nwl]  directional forward     scattering coefficient for diffuse  incidence
          cO.w[j,i]    = cO.sob[j]  * ρ_SW[j]   + cO.sof[j] *  τ_SW[j]; # [nl,nwl]  bidirectional scattering coefficent (directional-directional)
        end
    end
    cO.a    .= 1 .- cO.sigf;                  # [nl, nwl]     attenuation
end

# Compute optical properties of canopy
function computeCanopyGeomProps!(can::struct_canopy, angle::struct_angles, cO::struct_canopyOptProps)
    # Number of layers
    nl = nlayers
    # Leaf Area Index
    LAI = can.LAI
    dx = FT(1/nl)
    xl = collect(FT(0):FT(-1/nl):-1)
    tts = angle.tts
    tto = angle.tto
    psi = angle.psi
    lazitab = can.lazitab

    # Geometric quantities (ougoing direction first!)
    cto = cosd(tto)
    tanto	= tand(tto);
    cospsi	= cosd(psi);

    # Generate leaf angle distribution:
    if can.TypeLidf==1
        lidf,litab = dladgen(can.LIDFa,can.LIDFb);
    elseif can.TypeLidf==2
        # TBD Still need to check how to fix litab here:
        lidf,litab = campbell(can.LIDFa);
    end
    #println(lidf)

    # Precompute all of this before in a separate function?
    cos_ttlo    = cosd.(lazitab);  #   cos leaf azimuth angles
    cos_ttli    = cosd.(litab);    #   cosine of normal of upperside of leaf
    sin_ttli    = sind.(litab);    #   sine   of normal of upperside of leaf

    # Solar angle dependent ones:
    cts = cosd(tts)
    sin_tts  = sind(tts);
    ctscto = cts*cto;
    tants	= tand(tts);
    dso		= sqrt(tants*tants+tanto*tanto-2tants*tanto*cospsi);

    # angular distance, compensation of shadow length
    # Calculate geometric factors associated with extinction and scattering
    #Initialise sums
    cO.ks	= FT(0);
    cO.ko	= FT(0);
    cO.bf	= FT(0);
    cO.sob	= FT(0);
    cO.sof	= FT(0);

    #	Weighted sums over LIDF
	@inbounds for i=1:length(litab)
		# ttl = litab[i];	% leaf inclination discrete values
		ctl = cosd(litab[i]);
		#	SAIL volume scattering phase function gives interception and portions to be multiplied by rho and tau
		chi_s,chi_o,frho,ftau=volscatt(tts,tto,psi,litab[i]);
		#	Extinction coefficients
		ksli = chi_s./cts;
		koli = chi_o./cto;
        #	Area scattering coefficient fractions
		sobli	= frho*pi/ctscto;
		sofli	= ftau*pi/ctscto;
		bfli	= ctl*ctl;
		cO.ks	= cO.ks+ksli*lidf[i];
		cO.ko	= cO.ko+koli*lidf[i];
		cO.bf	= cO.bf+bfli*lidf[i];
		cO.sob	= cO.sob+sobli*lidf[i];
		cO.sof	= cO.sof+sofli*lidf[i];
    end
    #println(sob, " ", sof)

    #	Geometric factors to be used later with rho and tau
	cO.sdb	= (cO.ks+cO.bf)/2;
	cO.sdf	= (cO.ks-cO.bf)/2;
	cO.dob	= (cO.ko+cO.bf)/2;
	cO.dof	= (cO.ko-cO.bf)/2;
	cO.ddb	= (1 .+cO.bf)/2;
	cO.ddf	= (1 .-cO.bf)/2;

    # See eq 19 in vdT 2009
    Cs          = cos_ttli.*cts;             # [nli]     pag 305 modified by Joris
    Ss          = sin_ttli.*sin_tts;         # [nli]     pag 305 modified by Joris
    cos_deltas  = Cs*ones(FT,1,length(lazitab)) .+ Ss*cos_ttlo';   # [nli,nlazi]
    cO.fs          = abs.(cos_deltas./cts);         # [nli,nlazi] pag 305

    # 1.5 probabilities Ps, Po, Pso
    cO.Ps[:]          =   exp.(cO.ks*xl*LAI);     # [nl+1]  probability of viewing a leaf in solar dir
    cO.Po[:]          =   exp.(cO.ko*xl*LAI);     # [nl+1]  probability of viewing a leaf in observation dir

    cO.Ps[1:nl]    =   cO.Ps[1:nl] *(1 .-exp.(-cO.ks*LAI*dx))/(cO.ks*LAI*dx);   # Correct Ps/Po for finite dx
    cO.Po[1:nl]    =   cO.Po[1:nl] *(1 .-exp.(-cO.ko*LAI*dx))/(cO.ko*LAI*dx);   # Correct Ps/Po for finite dx

    #Pso: Probability of observing a sunlit leaf at depth x, see eq 31 in vdT 2009
    @inbounds for j=1:length(xl)
        cO.Pso[j] = quadgk(x -> Psofunction(cO.ko,cO.ks,LAI,can.hot,dso,x), xl[j]-dx,xl[j], rtol=1e-2)[1]/dx
    end

    #cO.Pso[cO.Pso.>cO.Po]= minimum([cO.Po[cO.Pso.>cO.Po] cO.Ps[cO.Pso.>cO.Po]],dims=2);    #takes care of rounding error
    #cO.Pso[cO.Pso.>cO.Ps]= minimum([cO.Po[cO.Pso.>cO.Ps] cO.Ps[cO.Pso.>cO.Ps]],dims=2);    #takes care of rounding error
end


# Had to make this like the Prosail F90 function, forwardDiff didn't work otherwise.
function calcJ1(x,m,k,LAI)
    del=(k-m)*LAI
    if abs(del)>1E-3;
        J1 = (exp(m*LAI*x)-exp(k*LAI*x))/(k-m);
    else
        J1 = -0.5*LAI*x*(exp(m*LAI*x)+exp(k*LAI*x))*(1.0-del^2/12.0);
    end
    return J1
end

function calcJ2(x,m,k,LAI)
    return (exp(k*LAI*x)-exp(-k*LAI)*exp(-m*LAI*(1+x)))/(k+m);
end


# APPENDIX IV function Pso from SCOPE v1.73
function Psofunction(K,k,LAI,q,dso,xl)
    if dso!=0.0
        alf         =   (dso/q) *2/(k+K);
        pso         =   exp((K+k)*LAI*xl + sqrt(K*k)*LAI/(alf  )*(1-exp(xl*(alf  ))));# [nl+1]  factor for correlation of Ps and Po
    else
        pso         =   exp((K+k)*LAI*xl - sqrt(K*k)*LAI*xl);# [nl+1]  factor for correlation of Ps and Po
    end
    return pso
end

# FROM SCOPE v1.73 APPENDIX II function volscat
"""********************************************************************************
!*	tts		= solar zenith
!*	tto		= viewing zenith
!*	psi		= azimuth
!*	ttl		= leaf inclination angle
!*	chi_s	= interception functions
!*	chi_o	= interception functions
!*	frho	= function to be multiplied by leaf reflectance rho
!*	ftau	= functions to be multiplied by leaf transmittance tau
!********************************************************************************

!	Compute volume scattering functions and interception coefficients
!	for given solar zenith, viewing zenith, azimuth and leaf inclination angle.

!	chi_s and chi_o are the interception functions.
!	frho and ftau are the functions to be multiplied by leaf reflectance rho and
!	leaf transmittance tau, respectively, in order to obtain the volume scattering
!	function.
"""
function volscatt(tts,tto,psi,ttli)
    #Volscatt version 2.
    #created by W. Verhoef
    #edited by Joris Timmermans to matlab nomenclature.
    # date: 11 February 2008
    #tts    [1]         Sun            zenith angle in degrees
    #tto    [1]         Observation    zenith angle in degrees
    #psi    [1]         Difference of  azimuth angle between solar and viewing position
    #ttli   [ttli]      leaf inclination array
    FT = typeof(tts);
    nli     = length(ttli);

    psi_rad         = deg2rad(psi);
    cos_psi         = cosd(psi);                #   cosine of relative azimuth angle
    cos_ttli        = cosd(ttli);               #   cosine of normal of upperside of leaf
    sin_ttli        = sind(ttli);               #   sine   of normal of upperside of leaf
    cos_tts         = cosd(tts);                #   cosine of sun zenith angle
    sin_tts         = sind(tts);                #   sine   of sun zenith angle
    cos_tto         = cosd(tto);                #   cosine of observer zenith angle
    sin_tto         = sind(tto);                #   sine   of observer zenith angle
    Cs              = cos_ttli*cos_tts;                 #   p305{1}
    Ss              = sin_ttli*sin_tts;                 #   p305{1}
    Co              = cos_ttli*cos_tto;                 #   p305{1}
    So              = sin_ttli*sin_tto;                 #   p305{1}
    As              = maximum([Ss,Cs]);
    Ao              = maximum([So,Co]);
    #println(-Cs./As, " ", Ss, " ", Cs)
    bts             = acos.(-Cs./As);                    #   p305{1}
    bto             = acos.(-Co./Ao);                    #   p305{2}
    chi_o           = 2/FT(pi)*((bto-FT(pi)/2).*Co + sin(bto).*So);
    chi_s           = 2/FT(pi)*((bts-FT(pi)/2).*Cs + sin(bts).*Ss);
    delta1          = abs(bts-bto);                     #   p308{1}
    delta2          = pi-abs(bts + bto - pi);           #   p308{1}

    Tot             = psi_rad + delta1 + delta2;        #   pag 130{1}

    bt1             = minimum([psi_rad,delta1]);
    bt3             = maximum([psi_rad,delta2]);
    bt2             = Tot - bt1 - bt3;

    T1              = 2Cs.*Co + Ss.*So.*cos_psi;
    T2              = sin(bt2).*(2As.*Ao + Ss.*So.*cos(bt1).*cos(bt3));

    Jmin            = (   bt2).*T1 - T2;
    Jplus           = (pi-bt2).*T1 + T2;

    frho            =  Jplus/(2pi^2);
    ftau            = -Jmin /(2pi^2);

    # pag.309 wl-> pag 135{1}
    frho            = maximum([FT(0),frho]);
    ftau            = maximum([FT(0),ftau]);
    #println(tts, " ",tto, " ",psi, " ",ttli)
    #println(chi_s, " ", chi_o, " ",frho, " ",ftau)
    return chi_s,chi_o,frho,ftau
end

function dladgen(a::Number,b::Number)

    freq = similar(litab)
    for i1=1:8
        t   = i1*10;
        freq[i1]=dcum(a,b,t);
    end
    for i2=9:12
        t   = 80.0+(i2-8)*2.;
        freq[i2]=dcum(a,b,t);
    end

    freq[13]=1;
    for i   = 13:-1:2
        freq[i]=freq[i]-freq[i-1];
    end
    return freq, litab
end

function dcum(a::Number,b::Number,t::Number)
    y = 0.0
    if a>=1
        f = 1-cosd(t);
    else
        epsi=1e-8;
        delx=1;
        x=2*deg2rad(t);
        p=x;
    	while (delx >= epsi)
            #println(delx)
            y = a*sin(x)+0.5*b*sin(2.0*x);
            dx=0.5*(y-x+p);
            x=x+dx;
            delx=abs(dx);
        end
    	f = (2.0*y+p)/pi;
    end
    return f
end

"""
From SCOPE v1.73:
********************************************************************************
*                          Campbell.f
*
*    Computation of the leaf angle distribution function value (freq)
*    Ellipsoidal distribution function caracterised by the average leaf
*    inclination angle in degree (ala)
*    Campbell 1986
*
********************************************************************************
 edit 2017 12 28: change sampling of angles to match with dladgen.m
"""
function campbell(ala)
    tx1=[10.,20.,30.,40.,50.,60.,70.,80.,82.,84.,86.,88.,90.];
    tx2=[0.,10.,20.,30.,40.,50.,60.,70.,80.,82.,84.,86.,88.];

    litab   = (tx2.+tx1)./2.0;
    n=length(litab);
    tl1     = deg2rad(tx1)
    tl2     = deg2rad(tx2)
    excent  = exp(-1.6184e-5*ala^3+2.1145e-3*ala^2-1.2390e-1*ala+3.2491);
    sum0    = 0;

    freq=zeros(n);
    for i=1:n
        x1  = excent./(sqrt(1 .+excent^2 .*tan(tl1(i)).^2));
        x2  = excent./(sqrt(1 .+excent^2 .*tan(tl2(i)).^2));
        if (excent==1)
            freq[i] = abs(cos(tl1(i))-cos(tl2(i)));
        else
            alpha  = excent./sqrt(abs(1-excent.^2));
            alpha2 = alpha.^2;
            x12 = x1.^2;
            x22 = x2.^2;
            if (excent>1)
                alpx1 = sqrt(alpha2(excent>1)+x12(excent>1));
                alpx2[excent>1] = sqrt(alpha2(excent>1)+x22(excent>1));
                dum   = x1*alpx1+alpha2*log(x1+alpx1);
                freq[i] = abs(dum-(x2.*alpx2+alpha2.*log(x2+alpx2)));
            else
                almx1 = sqrt(alpha2-x12);
                almx2 = sqrt(alpha2-x22);
                dum   = x1.*almx1+alpha2.*asin(x1./alpha);
                freq[i] = abs(dum-(x2.*almx2+alpha2.*asin(x2./alpha)));
            end
        end
    end
    sum0 = sum(freq,dims=2);
    freq0=freq./sum0;
return freq0,litab
end

"""
    calctav(alfa, nr)

    ***********************************************************************
    From calctav.m in PROSPECT-D
    ***********************************************************************
    Stern F. (1964), Transmission of isotropic radiation across an
    interface between two dielectrics, Appl. Opt., 3(1):111-113.
    Allen W.A. (1973), Transmission of isotropic light across a
    dielectric surface in two and three dimensions, J. Opt. Soc. Am.,
    63(6):664-666.
    ***********************************************************************

"""
function calctav(α,nr)
    tt = typeof(nr)
    rd  = tt(pi/180)
    n2  = nr^2
    np  = n2+1
    nm  = n2-1
    a   = (nr+1)*(nr+1)/2
    k   = -(n2-1)*(n2-1)/4
    sa  = sin(α*rd)
    if α!=90.0
        b1  = sqrt((sa^2-np/2).*(sa^2-np/2)+k)
    else
        b1 = 0
    end
    b2  = sa^2-np/2
    b   = b1-b2
    b3  = b^3
    a3  = a^3
    ts  = (k^2/(6*b3)+k/b-b/2)-(k^2/(6*a3)+k/a-a/2)
    tp1 = -2*n2*(b-a)/(np^2)
    tp2 = -2*n2*np*log(b/a)/(nm^2)
    tp3 = n2*(1/b-1/a)/2
    tp4 = 16*n2^2*(n2^2+1)*log((2*np*b-nm^2)/(2*np*a-nm^2))/(np^3*nm^2)
    tp5 = 16*n2^3*(1/(2*np*b-nm^2)-1/(2*np*a-nm^2))/(np^3)
    tp  = tp1+tp2+tp3+tp4+tp5
    tav = (ts+tp)/(2*sa^2)
    return tav
end

#function expint(x)
#    p = Poly([8.267661952366478e+00, -7.773807325735529e-01, -3.012432892762715e-01, -7.811863559248197e-02, -1.019573529845792e-02,-6.973790859534190e-04,-2.569498322115933e-05, -4.819538452140960e-07,  -3.602693626336023e-09])
#    polyv = polyval(p,real(x));
#    k = findall( abs(imag(x)) <= polyv );
#    -GSL.sf_expint_Ei(-x)
#end
# From Matlab!
p = Polynomials.Poly(convert(Array{FT},[8.267661952366478e+00, -7.773807325735529e-01, -3.012432892762715e-01, -7.811863559248197e-02, -1.019573529845792e-02,-6.973790859534190e-04,-2.569498322115933e-05, -4.819538452140960e-07,  -3.602693626336023e-09]))
const egamma=FT(0.57721566490153286061);
function expint(x::Number)

    polyv = Polynomials.polyval(p,real(x));
    if abs(imag(x)) <= polyv
        #initialization
        xk = x;
        yk = -egamma - log(xk);
        j = 1;
        pterm = xk;
        term = xk;
        while abs(term) > (eps(yk))
            yk = yk + term;
            j = j + 1;
            pterm = -xk.*pterm/j;
            term = pterm/j;
        end # end of the while loop
        y = yk;
    else
        n = FT(1);
        xk = x;
        am2 = FT(0)
        bm2 = FT(1)
        am1 = FT(1)
        bm1 = xk
        f = am1 / bm1;
        oldf = Inf;
        j = FT(2);

        while abs(f-oldf) > (100*eps(FT)*abs(f))
            alpha = n-1+(j/2); # note: beta= 1
            #calculate A(j), B(j), and f(j)
            a = am1 + alpha * am2;
            b = bm1 + alpha * bm2;

            # save new normalized variables for next pass through the loop
            #  note: normalization to avoid overflow or underflow
            am2 = am1 / b;
            bm2 = bm1 / b;
            am1 = a / b;
            bm1 = FT(1);

            f = am1;
            j = j+1;

            # calculate the coefficients for j odd
            alpha = (j-1)/2;
            beta = xk;
            a = beta * am1 + alpha * am2;
            b = beta * bm1 + alpha * bm2;
            am2 = am1 / b;
            bm2 = bm1 / b;
            am1 = a / b;
            bm1 = 1;
            oldf = f;
            f = am1;
            #println(typeof(xk), typeof(f), typeof(alpha))
            j = j+1;
        end

        y= exp(-xk) * f - 1im*FT(pi)*((real(xk)<0)&(imag(xk)==0));
    end
    return y
end
end
