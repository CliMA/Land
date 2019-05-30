module FluspectMod
#using GSL
using Polynomials
using Statistics
# Matlab reading
using MAT
# Numerical integration package (Simson rule)
using QuadGK

# Is this OK?
file_Opti = joinpath(dirname(pathof(FluspectMod)), "Optipar2017_ProspectD.mat")

const minwle  = 400.; # PAR range
const maxwle  = 700.;
const minwlf  = 650.; # SIF range
const maxwlf  = 850.;

# Doubling Adding layers
const ndub = 10

# Canopy Layers:
const nl = 60
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


# Enable downsampling spectral resolution to arbitrary grid (specify array with boundaries here)
# Don't like global arrays yet, need to make sure this is not creating performance issues.
function init(swl)
    global WL = swl
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

function fluspect(x::Vector; fqe::Number=0.0, Kab__=Kab, Kant__=Kant, KBrown__=KBrown, Kw__=Kw, Km__=Km, KcaV__=KcaV, KcaZ__=KcaZ)
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
kChlrel  = Cab*Kab./(Kall.*N.+eps());

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
kChl     =   kChlrel .* k;

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

    A11  = xf*Ih + Iv*xe';           A12 = (xf*xe').*(rf*Ih .+ Iv*re');
    A21  = 1 .+(xf*xe').*(1 .+rf*re');   A22 = (xf.*rf)*Ih+Iv*(xe.*re)';
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

function RTM_sail(x::Vector; TypeLidf=1)
    # State Vector X includes (in that order):
    #N,Cab,Car,Ant,Cbrown,Cw,Cm,Cx,LIDFa,LIDFb,lai,q,tts,tto,psi,rsoil
    LIDFa = x[9]
    LIDFb = x[10]
    LAI = x[11]
    q = x[12]
    tts = x[13]
    tto = x[14]
    psi = x[15]
    iLAI    = LAI/nl;               # [1] LAI of elementary layer (guess we can change that)

    #rsoil = x[16] # will be handled later

    LRT = fluspect(x[1:8], fqe=0.0)
    ρ=LRT[:,1]
    τ=LRT[:,2]
    # Geometric quantities (need to check allocation cost!)
    cts = cos(deg2rad(tts))
    cto = cos(deg2rad(tto))
    sin_tts  = sin(deg2rad(tts));             #           sin solar       angle
    ctscto = cts*cto;
    tants	= tan(deg2rad(tts));
    tanto	= tan(deg2rad(tto));
    cospsi	= cos(deg2rad(psi));


    dso		= sqrt(tants*tants+tanto*tanto-2.0*tants*tanto*cospsi);

    # Generate leaf angle distribution:
    if TypeLidf==1
        lidf,litab = dladgen(LIDFa,LIDFb);
    elseif TypeLidf==2
        lidf,litab = campbell(LIDFa);
    end
    #println(lidf)

    cos_ttlo    = cos.(deg2rad.(lazitab));  #   cos leaf azimuth angles
    cos_ttli    = cos.(deg2rad.(litab));    #   cosine of normal of upperside of leaf
    sin_ttli    = sin.(deg2rad.(litab));    #   sine   of normal of upperside of leaf

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
		ctl = cos(deg2rad(litab[i]));
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
    Cs          = cos_ttli.*cts;             # [nli]     pag 305 modified by Joris
    Ss          = sin_ttli.*sin_tts;         # [nli]     pag 305 modified by Joris

    cos_deltas  = Cs*ones(1,length(lazitab)) .+ Ss*cos_ttlo';   # [nli,nlazi]
    fs          = abs.(cos_deltas./cts);         # [nli,nlazi] pag 305

    # 1.5 probabilities Ps, Po, Pso
    Ps          =   exp.(ks*xl*LAI);                                              # [nl+1]  p154{1} probability of viewing a leaf in solar dir
    Po          =   exp.(ko*xl*LAI);                                              # [nl+1]  p154{1} probability of viewing a leaf in observation dir

    Ps[1:nl]    =   Ps[1:nl] *(1 .-exp.(-ks*LAI*dx))/(ks*LAI*dx);   # Correct Ps/Po for finite dx
    Po[1:nl]    =   Po[1:nl] *(1 .-exp.(-ko*LAI*dx))/(ko*LAI*dx);   # Correct Ps/Po for finite dx

    # Not yet sure what this is q is (canopy.hot), need to read up later:
    # canopy.hot  = canopy.leafwidth/canopy.hc;
    q           =   0.05;
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
    #	Interaction with the soil
    # Invent rsoil for now
    rsoil = zeros(length(rho_dd)).+0.1
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
    return [rso rsd rdd rdo]

end

function calcJ1(x,m,k,LAI)
    if abs(m-k)>1E-3;
        J1 = (exp(m*LAI*x)-exp(k*LAI*x))./(k-m);
    else
        J1 = -.5*(exp(m*LAI*x)+exp(k*LAI*x))*LAI.*x.*(1.0-1.0/12.0*(k-m).^2LAI^2x.^2);
    end
    return J1
end

function calcJ2(x,m,k,LAI)
    return (exp(k*LAI*x)-exp(-k*LAI)*exp(-m*LAI*(1+x)))./(k+m);
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
    cos_psi         = cos(deg2rad(psi));                #   cosine of relative azimuth angle
    cos_ttli        = cos(deg2rad(ttli));               #   cosine of normal of upperside of leaf
    sin_ttli        = sin(deg2rad(ttli));               #   sine   of normal of upperside of leaf
    cos_tts         = cos(deg2rad(tts));                #   cosine of sun zenith angle
    sin_tts         = sin(deg2rad(tts));                #   sine   of sun zenith angle
    cos_tto         = cos(deg2rad(tto));                #   cosine of observer zenith angle
    sin_tto         = sin(deg2rad(tto));                #   sine   of observer zenith angle
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
        f = 1-cos(deg2rad(t));
    else
        epsi=1e-8;
        delx=1;
        x=2*deg2rad(t);
        p=x;
    	while (delx >= epsi)
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
