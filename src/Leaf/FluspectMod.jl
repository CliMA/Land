module FluspectMod
#using GSL
using Polynomials
using Statistics
# Matlab reading
using MAT
# Numerical integration package (Simpson rule)
using QuadGK
using MathToolsMod

# Is this OK?
file_Opti = joinpath(dirname(pathof(FluspectMod)), "Optipar2017_ProspectD.mat")

# Struct for observation and solar angles
mutable struct angles{TT<:Number}
         tts::TT # Solar Zenith Angle in degrees
         tto::TT # Viewing Zenith Angle in degrees
         psi::TT # relatvie azimuth in degrees
end

# struct for canopy
mutable struct canopy{TT<:Number}
         nlayers::Integer # number of layers
         LAI::TT # Leaf Area Index
         hc::TT # relatie azimuth in degrees
end

const minwle  = 400.; # PAR range
const maxwle  = 700.;
const minwlf  = 650.; # SIF range
const maxwlf  = 850.;

# Doubling Adding layers
const ndub = 15

# Canopy Layers:
const nl = 20
const xl = collect(0:-1/nl:-1)
const lazitab = collect(5:10:355)
const dx = 1.0/nl

# Read in all optical data:
opti = matread(file_Opti)["optipar"]
nr_     =  convert(Array{Float32}, opti["nr"]);nr = nr_
Km_     =  convert(Array{Float32}, opti["Kdm"]);Km = Km_
Kab_    =  convert(Array{Float32}, opti["Kab"]);Kab = Kab_
Kant_   =  convert(Array{Float32}, opti["Kant"]);Kant =Kant_
Kcar_   =  convert(Array{Float32}, opti["Kca"]);Kcar= Kcar_
Kw_     =  convert(Array{Float32}, opti["Kw"]);Kw=Kw_
KBrown_ =  convert(Array{Float32}, opti["Ks"]);KBrown=KBrown_
phi_    =  convert(Array{Float32}, opti["phi"]);phi=phi_
KcaV_   =  convert(Array{Float32}, opti["KcaV"]);KcaV=KcaV_
KcaZ_   =  convert(Array{Float32}, opti["KcaZ"]);KcaZ =KcaZ_
lambda_ =  convert(Array{Float32}, opti["wl"]);lambda = lambda_
kChlrel = similar(lambda)

# Enable downsampling spectral resolution to arbitrary grid (specify array with boundaries here)
# Don't like global arrays yet, need to make sure this is not creating performance issues.
function init(swl)
    global WL = swl
    #Ipar    = find(wl>=minPAR & wl<=maxPAR);
    # Stupid...
    global nr = zeros(length(swl)-1)
    global Km = zeros(length(swl)-1)
    global Kab = zeros(length(swl)-1)
    global Kant = zeros(length(swl)-1)
    global Kcar = zeros(length(swl)-1)
    global Kw = zeros(length(swl)-1)
    global KBrown = zeros(length(swl)-1)
    global phi = zeros(length(swl)-1)
    global KcaV = zeros(length(swl)-1)
    global KcaZ = zeros(length(swl)-1)
    global lambda = zeros(length(swl)-1)
    global kChlrel = zeros(length(swl)-1)
    for i in 1:length(swl)-1
        wo = findall((lambda_.>=swl[i]).&(lambda_.<swl[i+1]) )
        #println(mean(nr_[wo]))
        nr[i]   =  mean(nr_[wo])
        Km[i]  =  mean(Km_[wo])
        Kab[i]  =  mean(Kab_[wo])
        Kant[i] =  mean(Kant_[wo])
        Kcar[i]  =  mean(Kcar_[wo])
        Kw[i]   =  mean(Kw_[wo])
        KBrown[i] = mean(KBrown_[wo])
        phi[i]  =  mean(phi_[wo])
        KcaV[i] = mean(KcaV_[wo])
        KcaZ[i] = mean(KcaZ_[wo])
        lambda[i] = mean(lambda_[wo])
    end
end

function fluspect!(x::Vector; fqe::Number=0.0, Kab__=Kab, Kant__=Kant, KBrown__=KBrown, Kw__=Kw, Km__=Km, KcaV__=KcaV, KcaZ__=KcaZ, kChlrel__ = kChlrel)
    #println(fqe)
# , Kab__=Kab, Kant__=Kant, KBrown__=KBrown, Kw__=Kw, Km__=Km, KcaV__=KcaV, KcaZ__=KcaZ)
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
N = x[1]
Cab = x[2]
Car = x[3]
Ant = x[4]
Brown = x[5]
Cw = x[6]
Cm = x[7]
Cx = x[8]
Kcaro = (1.0-Cx).* KcaV__ + Cx .* KcaZ__;

Kall    = (Cab*Kab__.+Car*Kcaro.+Ant*Kant__.+Brown*KBrown__.+Cw*Kw__.+Cm*Km__)./N
# Relative absorption by Chlorophyll only (drives SIF and GPP eventually)
kChlrel__  = Cab*Kab./(Kall.*N.+eps());

# Adding eps() here to keep it stable and NOT set to 1 manually when Kall=0 (ForwardDiff won't work otherwise)
tau = real.((1.0.-Kall).*exp.(-Kall) .+ Kall.^2.0.*expint.(Kall.+eps()))

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
talf    = calctav.(40.,nr)
ralf    = 1.0.-talf
t12     = calctav.(90.,nr)
r12     = 1.0.-t12
t21     = t12./(nr.^2)
r21     = 1.0.-t21

# top surface side
denom   = 1.0.-r21.*r21.*tau.^2
Ta      = talf.*tau.*t21./denom
Ra      = ralf.+r21.*tau.*Ta

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
D       = sqrt.((1.0.+r.+t).*(1.0.+r.-t).*(1.0.-r.+t).*(1.0.-r.-t))
rq      = r.^2
tq      = t.^2
a       = (1.0.+rq.-tq.+D)./(2r)
b       = (1.0.-rq.+tq.+D)./(2t)

bNm1    = b.^(N.-1);                  #
bN2     = bNm1.^2
a2      = a.^2
denom   = a2.*bN2.-1
Rsub    = a.*(bN2.-1)./denom
Tsub    = bNm1.*(a2.-1)./denom

# Case of zero absorption
j       = findall(r.+t .>= 1)
Tsub[j] = t[j]./(t[j]+(1.0.-t[j])*(N-1))
Rsub[j]	= 1.0.-Tsub[j]

# Reflectance & transmittance of the leaf: combine top layer with next N-1 layers
denom   = 1.0.-Rsub.*r
tran    = Ta.*Tsub./denom
refl    = Ra.+Ta.*Rsub.*t./denom

RT     = [refl tran]
if fqe ==0.0
    return RT
end
# FROM SCOPE notes:
# From here a new path is taken: The doubling method used to calculate
# fluoresence is now only applied to the part of the leaf where absorption
# takes place, that is, the part exclusive of the leaf-air interfaces. The
# reflectance (rho) and transmittance (tau) of this part of the leaf are
# now determined by "subtracting" the interfaces
# CF Note: All of the below takes about 10 times more time than the RT above. Need to rething speed and accuracy. (10nm is bringing it down a lot!)

Rb  = (refl.-ralf)./(talf.*t21+(refl.-ralf).*r21);  # Remove the top interface
Z   = tran.*(1.0.-Rb.*r21)./(talf.*t21);            # Derive Z from the transmittance

rho = (Rb.-r21.*Z.^2)./(1.0.-(r21.*Z).^2);      # Reflectance and transmittance
tau = (1.0.-Rb.*r21)./(1.0.-(r21.*Z).^2).*Z;    # of the leaf mesophyll layer
t   =   tau;
r   =   max.(rho,0.0);                       # Avoid negative r

# Derive Kubelka-Munk s and k
I_rt     =   findall((r.+t).<1);
D[I_rt]  =   sqrt.((1 .+ r[I_rt] .+ t[I_rt]) .* (1 .+ r[I_rt] .- t[I_rt]) .* (1 .- r[I_rt] .+ t[I_rt]) .*  (1 .- r[I_rt] .- t[I_rt]));
a[I_rt]  =   (1 .+ r[I_rt].^2 .- t[I_rt].^2 .+ D[I_rt]) ./ (2r[I_rt]);
b[I_rt]  =   (1 .- r[I_rt].^2 + t[I_rt].^2 .+ D[I_rt]) ./ (2t[I_rt]);
a[(r.+t).>=1] .=   1.0;
b[(r.+t).>=1] .=   1.0;

s        =   r./t;
I_a      =   findall((a.>1).&(a.!=Inf));
s[I_a]   =   2 .*a[I_a] ./ (a[I_a].^2 .- 1) .* log.(b[I_a]);

k        =   log.(b);
k[I_a]   =   (a[I_a].-1) ./ (a[I_a].+1) .* log.(b[I_a]);
kChl     =   kChlrel__ .* k;

# indices of wle and wlf within wlp
Iwle   = findall((lambda.>=minwle) .& (lambda.<=maxwle));
Iwlf   = findall((lambda.>=minwlf) .& (lambda.<=maxwlf));
wle    = lambda[Iwle];    # excitation wavelengths,  column
wlf    = lambda[Iwlf];    # fluorescence wavelengths,  column
epsi         = 2.0^(-ndub);

# initialisations
te     = 1 .-(k[Iwle].+s[Iwle]) .* epsi;
tf     = 1 .-(k[Iwlf].+s[Iwlf]) .* epsi;
re     = s[Iwle] .* epsi;
rf     = s[Iwlf] .* epsi;

sigmoid     = 1 ./(1 .+exp.(-wlf./10).*exp.(wle'./10));  # matrix computed as an outproduct
#println(size(sigmoid)," ", size(phi), " ", size(kChl)," ", size(Iwle), " ", size(Iwlf), " ", size(kChl[Iwle]))

Mf = Mb = fqe .* ((0.5*phi[Iwlf]).*epsi) .* kChl[Iwle]'.*sigmoid

Ih          = ones(1,length(te));     # row of ones
Iv          = ones(length(tf),1);     # column of ones

# Doubling routine
for i = 1:ndub
    xe = te./(1 .-re.*re);  ten = te.*xe;  ren = re.*(1 .+ten);
    xf = tf./(1 .-rf.*rf);  tfn = tf.*xf;  rfn = rf.*(1 .+tfn);

    A11  = xf*Ih + Iv*xe';
    A12 = (xf*xe').*(rf*Ih .+ Iv*re');
    A21  = 1 .+(xf*xe').*(1 .+rf*re');
    A22 = (xf.*rf)*Ih+Iv*(xe.*re)';
    #println(size(A11)," ", size(A12), " ", size(Mf)," ", size(Mb), " ")
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

Mb  = gn;
Mf  = fn;
return  RT,Mf,Mb
end

function RTM_sail(x::Vector; LIDFa=-0.35, LIDFb=-0.15,q=0.05, tts=30, tto=0, psi=90, wl=lambda[:,1], TypeLidf=1)
    # Index of PAR wavelengths
    iPAR = findall((wl.>minwle) .& (wl.<=maxwle))
    # State Vector X includes (in that order):
    #N,Cab,Car,Ant,Cbrown,Cw,Cm,Cx,lai,rsoil
    LAI = x[9]

    # Define soil as polynomial (depends on state vector size):
    pSoil =  Polynomials.Poly(x[10:end])
    rsoil = Polynomials.polyval(pSoil,wl.-mean(wl));

    iLAI    = LAI/nl;               # [1] LAI of elementary layer (guess we can change that)

    # Size of wavelength array:
    nwl = length(wl)

    # Call Fluspect for Leaf optical properties (can be done outside later if needed, need to test timing)
    LRT = fluspect!(x[1:8], fqe=0.0)
    ρ=LRT[:,1]
    τ=LRT[:,2]
    leaf_emiss = 1.0.-ρ.-τ
    # Geometric quantities (need to check allocation cost!)
    cts = cosd(tts)
    cto = cosd(tto)
    sin_tts  = sind(tts);             #           sin solar       angle
    ctscto = cts*cto;
    tants	= tand(tts);
    tanto	= tand(tto);
    cospsi	= cosd(psi);


    dso		= sqrt(tants*tants+tanto*tanto-2.0*tants*tanto*cospsi);

    # Generate leaf angle distribution:
    if TypeLidf==1
        lidf,litab = dladgen(LIDFa,LIDFb);
    elseif TypeLidf==2
        lidf,litab = campbell(LIDFa);
    end
    #println(lidf)

    cos_ttlo    = cosd.(lazitab);  #   cos leaf azimuth angles
    cos_ttli    = cosd.(litab);    #   cosine of normal of upperside of leaf
    sin_ttli    = sind.(litab);    #   sine   of normal of upperside of leaf

    # angular distance, compensation of shadow length
    # Calculate geometric factors associated with extinction and scattering
    #Initialise sums
    ks	= 0.0;
    ko	= 0.0;
    bf	= 0.0;
    sob	= 0.0;
    sof	= 0.0;

    #	Weighted sums over LIDF
	for i=1:length(litab)
		# ttl = litab[i];	% leaf inclination discrete values
		ctl = cosd(litab[i]);
		#	SAIL volume scattering phase function gives interception and portions to be
		#	multiplied by rho and tau
		chi_s,chi_o,frho,ftau=volscatt(tts,tto,psi,litab[i]);

		#********************************************************************************
		#*                   SUITS SYSTEM COEFFICIENTS
		#*
		#*	ks  : Extinction coefficient for direct solar flux
		#*	ko  : Extinction coefficient for direct observed flux
		#*	att : Attenuation coefficient for diffuse flux
		#*	sigb : Backscattering coefficient of the diffuse downward flux
		#*	sigf : Forwardscattering coefficient of the diffuse upward flux
		#*	sf  : Scattering coefficient of the direct solar flux for downward diffuse flux
		#*	sb  : Scattering coefficient of the direct solar flux for upward diffuse flux
		#*	vf   : Scattering coefficient of upward diffuse flux in the observed direction
		#*	vb   : Scattering coefficient of downward diffuse flux in the observed direction
		#*	w   : Bidirectional scattering coefficient
		#********************************************************************************

		#	Extinction coefficients
		ksli = chi_s./cts;
		koli = chi_o./cto;

		#	Area scattering coefficient fractions
		sobli	= frho*pi/ctscto;
		sofli	= ftau*pi/ctscto;
		bfli	= ctl*ctl;
		ks	= ks+ksli*lidf[i];
		ko	= ko+koli*lidf[i];
		bf	= bf+bfli*lidf[i];
		sob	= sob+sobli*lidf[i];
		sof	= sof+sofli*lidf[i];
    end
    #println(sob, " ", sof)
    #	Geometric factors to be used later with rho and tau
	sdb	= 0.5*(ks+bf);
	sdf	= 0.5*(ks-bf);
	dob	= 0.5*(ko+bf);
	dof	= 0.5*(ko-bf);
	ddb	= 0.5*(1 .+bf);
	ddf	= 0.5*(1 .-bf);
    # Skipped SCOPE lines 186-213 here (catch up later)
    # 1.4 solar irradiance factor for all leaf orientations
    # See eq 19 in vdT 2009
    Cs          = cos_ttli.*cts;             # [nli]     pag 305 modified by Joris
    Ss          = sin_ttli.*sin_tts;         # [nli]     pag 305 modified by Joris

    cos_deltas  = Cs*ones(1,length(lazitab)) .+ Ss*cos_ttlo';   # [nli,nlazi]
    fs          = abs.(cos_deltas./cts);         # [nli,nlazi] pag 305

    # 1.5 probabilities Ps, Po, Pso
    Ps          =   exp.(ks*xl*LAI);     # [nl+1]  p154{1} probability of viewing a leaf in solar dir
    Po          =   exp.(ko*xl*LAI);     # [nl+1]  p154{1} probability of viewing a leaf in observation dir

    Ps[1:nl]    =   Ps[1:nl] *(1 .-exp.(-ks*LAI*dx))/(ks*LAI*dx);   # Correct Ps/Po for finite dx
    Po[1:nl]    =   Po[1:nl] *(1 .-exp.(-ko*LAI*dx))/(ko*LAI*dx);   # Correct Ps/Po for finite dx



    #Pso: Probability of observing a sunlit leaf at depth x, see eq 31 in vdT 2009
    Pso         =   similar(Po);
    for j=1:length(xl)
        #println(size(a), " ", size(Pso), " ", size(Po))
        Pso[j] = quadgk(x -> Psofunction(ko,ks,LAI,q,dso,x), xl[j]-dx,xl[j], rtol=1e-2)[1]/dx
        #Pso[j,:]=   quad(@(y)Psofunction(K,k,LAI,q,dso,y),xl(j)-dx,xl(j))/dx; %#ok<FREMO>
    end

    Pso[Pso.>Po]= minimum([Po[Pso.>Po] Ps[Pso.>Po]],dims=2);    #takes care of rounding error
    Pso[Pso.>Ps]= minimum([Po[Pso.>Ps] Ps[Pso.>Ps]],dims=2);    #takes care of rounding error


    # All with length of wavelengths:
    sigb = ddb.*ρ.+ddf.*τ;
	sigf = ddf.*ρ.+ddb.*τ;
    sb   = sdb*ρ .+ sdf*τ;            # [nwl]     sb,     p305{1} diffuse     backscatter scattering coefficient for specular incidence
    sf   = sdf*ρ .+ sdb*τ;            # [nwl]     sf,     p305{1} diffuse     forward     scattering coefficient for specular incidence
    vb   = dob*ρ .+ dof*τ;            # [nwl]     vb,     p305{1} directional backscatter scattering coefficient for diffuse  incidence
    vf   = dof*ρ .+ dob*τ;            # [nwl]     vf,     p305{1} directional forward     scattering coefficient for diffuse  incidence
    w    = sob*ρ .+ sof*τ;            # [nwl]     w,      p309{1} bidirectional scattering coefficent (directional-directional)
    a    = 1 .-sigf;                  # [nwl]     attenuation
    m    = sqrt.(a.^2 .-sigb.^2);     # [nwl]
    rinf = (a.-m)./sigb;              # [nwl]
    rinf2= rinf.*rinf;                # [nwl]
    #println(minimum(m), " ", min(ks))
    # direct solar radiation
    J1k        = calcJ1.(-1, m,ks,LAI);          # [nwl]
    J2k        = calcJ2.( 0, m,ks,LAI);          # [nwl]
    J1K        = calcJ1.(-1, m,ko,LAI);          # [nwl]   % added for calculation of rdo
    J2K        = calcJ2.( 0, m,ko,LAI);          # [nwl]   % added for calculation of rdo
    e1          = exp.(-m.*LAI);                 # [nwl]
    e2          = e1.^2;                         # [nwl]
    re          = rinf.*e1;                      # [nwl]

    denom       = 1 .-rinf2.*e2;
    s1          = sf .+rinf.*sb;
    s2          = sf.*rinf+sb;
    v1          = vf.+rinf.*vb;
    v2          = vf.*rinf.+vb;

    Pss         = s1.*J1k;          # [nwl]
    Qss         = s2.*J2k;          # [nwl]

    Poo         = v1.*J1K;          # (nwl)   % added for calculation of rdo
    Qoo         = v2.*J2K;          # [nwl]   % added for calculation of rdo

    tau_ss      = exp(-ks*LAI);                  # [1]
    tau_oo      = exp(-ko*LAI);

    Z           = (1 - tau_ss * tau_oo)/(ks + ko);  # needed for analytic rso

    tau_dd      = (1 .-rinf2).*e1 ./denom;        # [nwl]
    rho_dd      = rinf.*(1 .-e2)  ./denom;        # [nwl]
    tau_sd      = (Pss.-re.*Qss) ./denom;         # [nwl]
    tau_do      = (Poo.-re.*Qoo) ./denom;         # [nwl]
    rho_sd      = (Qss.-re.*Pss) ./denom;         # [nwl]
    rho_do      = (Qoo.-re.*Poo) ./denom;         # (nwl)

    T1          = v2.*s1.*(Z.-J1k*tau_oo)./(ko.+m).+v1.*s2.*(Z.-J1K*tau_ss)./(ks.+m);
    T2          = -(Qoo.*rho_sd+Poo.*tau_sd).*rinf;
    rho_sod     = (T1+T2)./(1 .-rinf2);

    #	Bidirectional reflectance
    #	Single scattering contribution
    rho_sos = w.*iLAI.*sum(Pso[1:nl]);
    #	Total canopy contribution
    rho_so=rho_sos.+rho_sod;
    #println(rho_so[100:120])

    dn=1 .-rsoil.*rho_dd;
    # Total canopy contribution
    rso         = rho_so .+ rsoil .* Pso[nl+1] .+ ((tau_sd.+tau_ss*rsoil.*rho_dd).*tau_oo.+(tau_sd.+tau_ss).*tau_do).*rsoil./denom;
    # SAIL analytical reflectances
    # rsd: directional-hemispherical reflectance factor for solar incident flux
    rsd     = rho_sd .+ (tau_ss .+ tau_sd).*rsoil.*tau_dd./denom;
    # rdd: bi-hemispherical reflectance factor
    rdd     = rho_dd .+ tau_dd.*rsoil.*tau_dd./denom;
    # rdo: bi-directional reflectance factor
    rdo     = rho_do .+ (tau_oo .+ tau_do).*rsoil.*tau_dd./denom;
    #return [rso rsd rdd rdo]

    # Dummy code here to track direct and diffuse light first, need to separate this into another function later
    Esun_ = zeros(nwl).+200
    Esky_ = zeros(nwl).+100
    Emin_ = zeros(nl+1,nwl)
    Eplu_ = zeros(nl+1,nwl)

    Eplu_1      = rsoil.*((tau_ss.+tau_sd).*Esun_.+tau_dd.*Esky_)./denom;
    Eplu0       = rho_sd.*Esun_ .+ rho_dd.*Esky_ .+ tau_dd.*Eplu_1;
    Emin_1      = tau_sd.*Esun_ .+ tau_dd.*Esky_ .+ rho_dd.*Eplu_1;
    delta1      = Esky_  .- rinf.*Eplu0;
    delta2      = Eplu_1 .- rinf.*Emin_1;
    #println("Eplu1 ", Eplu_1[1:10])
    #println("Eplu0 ", Eplu0[1:10])
    #println("Emin_1 ", Emin_1[1:10])
    #println("delta1 ", delta1[1:10])
    #println("delta1 ", delta2[1:10])
    # calculation of the fluxes in the canopy (this seems slow!)
    t1 = sf.+rinf.*sb
    t2 = sb+rinf.*sf

    # The order here mattered, now doing the loop over nl, not nwl! (faster)
    # This loop is time consuming, probably don't need high spectral resolution for the NIR part here, so it can be shortened a lot (say 100nm steps from 700-2500nm?)
    # We just need it for net energy balance and PAR here.
    @inbounds for i = 1:nl+1
        J1kx  = calcJ1.(xl[i],m,ks,LAI);    #           [nl]
        J2kx  = calcJ2.(xl[i],m,ks,LAI);    #           [nl]
        F1    = Esun_.*J1kx.*t1 .+ delta1.*exp.( m.*LAI.*xl[i]);      #[nl]
        F2    = Esun_.*J2kx.*t2 .+ delta2.*exp.(-m.*LAI.*(xl[i].+1));  #[nl]
        Emin_[i,:]  = (F1.+rinf.*F2)./(1 .-rinf2);#        [nl,nwl]
        Eplu_[i,:]  = (F2.+rinf.*F1)./(1 .-rinf2);#        [nl,nwl]
    end

    #println(length(wl[iPAR]))
    Psun        = 0.001 * simpson_nu(wl[iPAR], e2phot(wl[iPAR]*1E-9,Esun_[iPAR]));
    #Psky        = 0.001 * helpers.Sint(e2phot(wlPAR*1E-9,Esky_(Ipar)),wlPAR);
    Asun        = 0.001 * simpson_nu(Esun_.*leaf_emiss,wl);                         # Total absorbed solar radiation
    Pnsun       = 0.001 * simpson_nu(wl[iPAR],e2phot(wl[iPAR]*1E-9,Esun_[iPAR].*leaf_emiss[iPAR]));  # Absorbed solar radiation  in PAR range in moles m-2 s-1
    Rnsun_PAR   = 0.001 * simpson_nu(wl[iPAR], Esun_[iPAR].*leaf_emiss[iPAR]);
    Pnsun_Cab   = 0.001 * simpson_nu(wl[iPAR], e2phot(wl[iPAR]*1E-9,kChlrel[iPAR].*Esun_[iPAR].*leaf_emiss[iPAR]));
    #println("Psun ", Psun)

    # 3. outgoing fluxes, hemispherical and in viewing direction, spectrum
    # hemispherical, spectral (Eout_)
    Eout_   = Eplu_[1,:];                  #           [nwl]

    # in viewing direction, spectral
    piLoc_      = (vb.*(Emin_[1:nl,:]'*Po[1:nl]) + vf.*(Eplu_[1:nl,:]'*Po[1:nl]) + w.*Esun_*sum(Pso[1:nl]))*iLAI;
    #println(size(rsoil)," ", size(Emin_[nl+1,:])," ", Po[nl+1], " ",Pso[nl+1])
    piLos_      = rsoil.*Emin_[nl+1,:].*Po[nl+1] + rsoil.*Esun_.*Pso[nl+1];
    piLo_       = piLoc_ + piLos_;              #           [nwl]
    Lo_         = piLo_/pi;

    # 4. net fluxes, spectral and total, and incoming fluxes
    # incident PAR at the top of canopy, spectral and spectrally integrated
    P_          = e2phot(wl[iPAR]*1E-9,(Esun_[iPAR]+Esky_[iPAR]));
    P           = .001 * simpson_nu(wl[iPAR], P_);


    # total direct radiation (incident and net) per leaf area (W m-2 leaf)
    Pdir        = fs * Psun;                        # [13 x 36]   incident
    Rndir       = fs * Asun;                        # [13 x 36]   net
    Pndir       = fs * Pnsun;                       # [13 x 36]   net PAR
    Pndir_Cab   = fs * Pnsun_Cab;                   # [13 x 36]   net PAR Cab
    Rndir_PAR   = fs * Rnsun_PAR;                   # [13 x 36]   net PAR energy units

    # Allocate memory (should pre-allocate in update and provide in structure!)
    Pdif = zeros(nl);
    Rndif= zeros(nl);
    Pndif= zeros(nl);
    Pndif_Cab= zeros(nl);
    Rndif_PAR= zeros(nl);
    Rndif_ = zeros(nl,nwl);
    Pndif_ = zeros(nl,length(iPAR));
    Pndif_Cab_ = zeros(nl,length(iPAR));
    Rndif_PAR_ = zeros(nl,length(iPAR));

    # canopy layers, diffuse radiation
    @inbounds for j = 1:nl
        # diffuse incident radiation for the present layer 'j' (mW m-2 um-1)
        E_         = .5*(Emin_[j,:] + Emin_[j+1,:]+ Eplu_[j,:]+ Eplu_[j+1,:]);

        # incident PAR flux, integrated over all wavelengths (moles m-2 s-1)
        Pdif[j]    = .001 * simpson_nu(wl[iPAR], e2phot(wl[iPAR]*1E-9,E_[iPAR]));  # [nl] , including conversion mW >> W

        # net radiation (W m-2) and net PAR (moles m-2 s-1), integrated over all wavelengths
        Rndif[j]           = .001 * simpson_nu(wl, Rndif_[j,:]);              # [nl]  Full spectrum net diffuse flux
        Pndif[j]           =        simpson_nu(wl[iPAR],Pndif_[j,iPAR]);      # [nl]  Absorbed PAR
        Pndif_Cab[j]       =        simpson_nu(wl[iPAR],Pndif_Cab_[j,iPAR]);  # [nl]  Absorbed PAR by Cab integrated
        Rndif_PAR[j]       = .001 * simpson_nu(wl[iPAR],Rndif_PAR_[j,iPAR]);  # [nl]  Absorbed PAR by Cab integrated

        # net radiation (mW m-2 um-1) and net PAR (moles m-2 s-1 um-1), per wavelength
        Rndif_[j,:]         = E_.*leaf_emiss;                                                   # [nl,nwl]  Net (absorbed) radiation by leaves
        Pndif_[j,:]         = .001 *(e2phot(wl[iPAR]*1E-9, Rndif_[j,iPAR]));             # [nl,nwlPAR]  Net (absorbed) as PAR photons
        Pndif_Cab_[j,:]     = .001 *(e2phot(wl[iPAR]*1E-9, kChlrel[iPAR].*Rndif_[j,iPAR]));  # [nl,nwl]  Net (absorbed) as PAR photons by Cab
        Rndif_PAR_[j,:]     = Rndif_[j,iPAR];  # [nl,nwlPAR]  Net (absorbed) as PAR energy
    end
    #Esunto  = 0.001 * simpson_nu(wl, Esun_) #Calculate optical sun fluxes (by Integration), including conversion mW >> W
    #Eskyto  = 0.001 * simpson_nu(wl, Esky_) #Calculate optical sun fluxes (by Integration)
    #println("Esunto ", Esunto)
    #println("Eskyto ", Eskyto)
    return Eplu_
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
    chi_o           = 2/pi*((bto-pi/2).*Co + sin(bto).*So);
    chi_s           = 2/pi*((bts-pi/2).*Cs + sin(bts).*Ss);
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
    frho            = maximum([0.0,frho]);
    ftau            = maximum([0.0,ftau]);
    #println(tts, " ",tto, " ",psi, " ",ttli)
    #println(chi_s, " ", chi_o, " ",frho, " ",ftau)
    return chi_s,chi_o,frho,ftau
end

function dladgen(a::Number,b::Number)
    litab=[5.,15.,25.,35.,45.,55.,65.,75.,81.,83.,85.,87.,89.];
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
    return freq,litab
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
    rd  = pi/180
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
    ts  = (k^2.0/(6*b3)+k/b-b/2)-(k^2.0/(6*a3)+k/a-a/2)
    tp1 = -2*n2*(b-a)/(np^2)
    tp2 = -2*n2*np*log(b/a)/(nm^2)
    tp3 = n2*(1.0/b-1.0/a)/2
    tp4 = 16*n2^2.0*(n2^2+1)*log((2*np*b-nm^2)/(2*np*a-nm^2))/(np^3.0*nm^2)
    tp5 = 16*n2^3.0*(1.0/(2*np*b-nm^2)-1.0/(2*np*a-nm^2))/(np^3)
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
function expint(x)
    p = Polynomials.Poly([8.267661952366478e+00, -7.773807325735529e-01, -3.012432892762715e-01, -7.811863559248197e-02, -1.019573529845792e-02,-6.973790859534190e-04,-2.569498322115933e-05, -4.819538452140960e-07,  -3.602693626336023e-09])
    polyv = Polynomials.polyval(p,real(x));
    if abs(imag(x)) <= polyv
        #initialization
        egamma=0.57721566490153286061;
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
        n = 1.0;
        xk = x;
        am2 = 0.0
        bm2 = 1.0
        am1 = 1.0
        bm1 = xk
        f = am1 / bm1;
        oldf = Inf;
        j = 2;
        while abs(f-oldf) > (100*eps()*abs(f))
            alpha = n-1+(j/2); # note: beta= 1
            #calculate A(j), B(j), and f(j)
            a = am1 + alpha * am2;
            b = bm1 + alpha * bm2;

            # save new normalized variables for next pass through the loop
            #  note: normalization to avoid overflow or underflow
            am2 = am1 / b;
            bm2 = bm1 / b;
            am1 = a / b;
            bm1 = 1.0;

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
            j = j+1;
        end
        y= exp(-xk) * f - 1im*pi*((real(xk)<0)&(imag(xk)==0));
    end
    return y
end


end
