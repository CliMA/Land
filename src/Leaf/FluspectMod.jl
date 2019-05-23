module FluspectMod
#using GSL
using Polynomials
using Statistics
# Matlab reading
using MAT

# This needs to be changed, more flexible instead of hard-coded:
const file_Opti = "/home/cfranken/code/gitHub/LSM-SPAM/src/Leaf/Optipar2017_ProspectD.mat"
const minwle  = 400.; # PAR range
const maxwle  = 700.;
const minwlf  = 650.; # SIF range
const maxwlf  = 850.;

# Doubling Adding layers
const ndub = 10

# Read in all optical data:
opti = matread(file_Opti)["optipar"]
nr_     =  opti["nr"];nr = nr_
Km_     =  opti["Kdm"];Km = Km_
Kab_    =  opti["Kab"];Kab = Kab_
Kant_   =  opti["Kant"];Kant =Kant_
Kcar_   =  opti["Kca"];Kcar= Kcar_
Kw_     =  opti["Kw"];Kw=Kw_
KBrown_ =  opti["Ks"];KBrown=KBrown_
phi_    =  opti["phi"];phi=phi_
KcaV_   =  opti["KcaV"];KcaV=KcaV_
KcaZ_   =  opti["KcaZ"];KcaZ =KcaZ_
lambda_ =  opti["wl"];lambda = lambda_


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

function fluspect(x::Vector;fqe::Number=0.01)
    println(fqe)
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
Kcaro = (1.0-Cx).* KcaV + Cx .* KcaZ;

Kall    = (Cab*Kab.+Car*Kcaro.+Ant*Kant.+Brown*KBrown.+Cw*Kw.+Cm*Km)./N
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
return  RT,Mf Mb
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
