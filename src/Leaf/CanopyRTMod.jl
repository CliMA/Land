module CanopyRTMod
#using GSL
using Polynomials
using Statistics
using ..PhysCon
using Parameters
# Matlab reading
using MAT
# Numerical integration package (Simpson rule)
using QuadGK

using ..Leaf


export optipar,
       angles,
       leafbio,
       leaf,
       optis,
       canopy,
       canRad,
       canOpt,
       FT,
       sunRad

# Floating point precision can be changed here
const FT = Float32


include("photo_structs.jl")

# Leaf inclination distribution
const litab      = FT[5.,15.,25.,35.,45.,55.,65.,75.,81.,83.,85.,87.,89.];
const litab_bnd  = FT[[0.,10.,20.,30.,40.,50.,60.,70.,80.,82.,84.,86.,88.] [10.,20.,30.,40.,50.,60.,70.,80.,82.,84.,86.,88., 90.] ];

const minwlPAR  = FT(400.); # PAR range
const maxwlPAR  = FT(700.);

const minwle  = FT(400.); # SIF excitation range
const maxwle  = FT(750.);

const minwlf  = FT(640.); # SIF emission range
const maxwlf  = FT(850.);

const fac     = FT(0.001) #(unit conversion factor mW >> W)

# Doubling Adding layers
const ndub = 10

# Load some standard settings here (can be done differently later but not bad to be fixed here)
# Wavelength boundaries:
#const swl = collect(FT(400):FT(5):2401)
const swl = [collect(FT(400):FT(10):650); collect(FT(655):FT(5):770); collect(FT(780):FT(25):2401)]
const dwl = diff(swl)
#println(dwl)

const optis = optipar{FT, length(swl)-1}(loadOpti(swl)...);

# Set wavelength here (will stay constant!!)
const wl   = optis.lambda
const nwl  = length(wl)
const Iwle = findall((wl.>=minwle) .& (wl.<=maxwle));
const Iwlf = findall((wl.>=minwlf) .& (wl.<=maxwlf));
const nWlF = length(Iwlf)
const iPAR = findall((wl.>=minwlPAR) .& (wl.<=maxwlPAR))
const wle  = wl[Iwle];    # excitation wavelengths,  column
const wlf  = wl[Iwlf];    # fluorescence wavelengths,  column

# Later on, we can implement an array of these structs for global parallel calculations!
# Load sun (will be coming from CLIMA in future, just file for now)
sunRad = incomingRadiation{FT}(loadSun(swl)...)

# Lots of structures predefine here, in principle just placeholders, routines can be run with any!
# Soil structure (can adapt later, just using a plain albedo here)
# WL, Shortwave_albedo, LE_albedo, Tskin:
const soil   = struct_soil{FT}(wl, 0.4*ones(FT, length(wl)), [0.1],290.0)
# Leaf structure (need to define the dimensions as in the other structures soon, tbd)
const leaf   = leafbio{FT,nwl, length(wle), length(wlf),length(wle)*length(wlf)}()
# Observation angles
const angles = struct_angles{FT}()
# Canopy Structure
const canopy = struct_canopy{FT}()
# Canopy radiation:
const canRad = struct_canopyRadiation{FT,nwl,nWlF,length(litab), length(canopy.lazitab), canopy.nlayers}()
# Canopy Optical Properties:
#const canOpt = struct_canopyOptProps{FT,nwl,canopy.nlayers, length(canopy.lazitab), length(litab)}()
canOpt = create_canopyOpt(FType=FT,nWL=nwl,nLayers=canopy.nlayers, nAzi=length(canopy.lazitab), nIncl=length(litab))

const nlayers = canopy.nlayers

function RTM_SW!(can::struct_canopy, cO::struct_canopyOptProps, cR::struct_canopyRadiation, sun::incomingRadiation, so::struct_soil)
    # Unpack variables from can structure
    @unpack Ω,nlayers, LAI = can


    iLAI    = LAI*Ω/nlayers;
    # Define soil as polynomial (depends on state vector size), still TBD in the structure mode now.
    rsoil = so.albedo_SW;
    soil_emissivity = 1 .-rsoil # emissivity of soil (goes into structure later)

    # 2.2. convert scattering and extinction coefficients into thin layer reflectances and transmittances
    # Eq. 17 in mSCOPE paper (changed here to compute real transmission)
    τ_ss = FT(exp(-cO.ks*iLAI));
    # Here, transmission looks ok, as it has to be sigf for iLAI=1!
    τ_dd = 1 .-cO.a*iLAI;
    τ_sd = cO.sf*iLAI;
    ρ_sd = cO.sb*iLAI;
    ρ_dd = cO.sigb*iLAI;

    # 2.3. reflectance calculation
    # Eq. 18 in mSCOPE paper
    # Just a scalar here but better to write it out:
    Xss  = τ_ss;

    # Soil reflectance boundary condition (same for diffuse and direct)
    cO.R_sd[:,end] = rsoil
    cO.R_dd[:,end] = rsoil

    # Eq. 18 in mSCOPE paper
    # from bottom to top:
    @inbounds for j=nlayers:-1:1
        dnorm       = 1 .-ρ_dd[:,j].*cO.R_dd[:,j+1];
        cO.Xsd[:,j]    = (τ_sd[:,j]+Xss*cO.R_sd[:,j+1].*ρ_dd[:,j])./dnorm;
        cO.Xdd[:,j]    = τ_dd[:,j]./dnorm;
        cO.R_sd[:,j]   = ρ_sd[:,j]+τ_dd[:,j].*(cO.R_sd[:,j+1]*Xss+cO.R_dd[:,j+1].*cO.Xsd[:,j]);
        cO.R_dd[:,j]   = ρ_dd[:,j]+τ_dd[:,j].*cO.R_dd[:,j+1].*cO.Xdd[:,j];
    end

    # 3.2 flux profile calculation
    # Eq. 19 in mSCOPE paper
    # Boundary condition at top: Incoming solar radiation
    #sun.E_diffuse = sun.E_diffuse.+0.2# Add some here for checking, is almost 0 otherwise!
    cO.Es_[:,1]    = sun.E_direct;
    cR.E_down[:,1] = sun.E_diffuse;

    @inbounds for j=1:nlayers # from top to bottom
        cR.netSW_sunlit[:,j] = cO.Es_[:,j].*(1 .-(τ_ss.+τ_sd[:,j]+ρ_sd[:,j]));
        cO.Es_[:,j+1]    =  Xss.*cO.Es_[:,j];
        cR.E_down[:,j+1] =  cO.Xsd[:,j].*cO.Es_[:,j]+cO.Xdd[:,j].*cR.E_down[:,j];
        cR.E_up[:,j]     =  cO.R_sd[:,j].*cO.Es_[:,j]+cO.R_dd[:,j].*cR.E_down[:,j];
    end
    # Boundary condition at the bottom, soil reflectance (Lambertian here)
    cR.E_up[:,end] = cO.R_sd[:,end].*cO.Es_[:,end]+cO.R_dd[:,end].*cR.E_down[:,end];

    # Hemispheric total outgoing (needs to go into stucture later):
    cR.Eout[:] = cR.E_up[:,1]

    # compute net diffuse radiation per layer:
    @inbounds for j=1:nlayers
        E_ = cR.E_down[:,j] +  cR.E_up[:,j+1];
        cR.netSW_shade[:,j]= E_.*(1 .-(τ_dd[:,j]+ρ_dd[:,j]))
        # Add diffuse radiation to direct radiation as well:
        #cR.netSW_sunlit[:,j] += cR.netSW_shade[:,j]
    end

    # outgoing in viewing direction
    # From Canopy
    piLoc_  = iLAI*sum(cO.vb.*cO.Po[1:nlayers]'.*cR.E_down[:,1:nlayers]+cO.vf.*cO.Po[1:nlayers]'.*cR.E_up[:,1:nlayers]+cO.w.*cO.Pso[1:nlayers]'.*sun.E_direct,dims=2)[:,1];
    # From Soil
    piLos_  = cR.E_up[:,end]*cO.Po[end];
    piLo_   = piLoc_ + piLos_;
    cR.Lo   = piLo_/pi;
    # Save albedos (hemispheric direct and diffuse as well as directional (obs))
    cR.alb_obs = piLo_./(sun.E_direct+sun.E_diffuse); # rso and rdo are not computed separately
    cR.alb_direct  = cO.R_sd[:,1]
    cR.alb_diffuse = cO.R_dd[:,1]

end;

function RTM_diffuseS(τ_dd::Array, ρ_dd::Array,S⁻::Array, S⁺::Array, boundary_top::Array, boundary_bottom::Array, rsoil::Array)
    FT = eltype(τ_dd)
    #Get dimensions (1st is wavelength, 2nd is layers), for Stefab Boltzmann, just one effective wavelength
    nwl,nl = size(τ_dd);
    Xdd = similar(τ_dd);
    Rdd = zeros(FT,nwl,nl+1)
    # Create Y,U matrices (see mSCOPE paper, eqs 30-39)
    Y  = similar(τ_dd);
    U  = similar(Rdd);
    E⁻ = similar(Rdd);
    E⁺ = similar(Rdd);
    net_diffuse = similar(τ_dd);

    Rdd[:,end] = rsoil

    # Source at TOC (0 for SIF, downwelling LW from atmosphere for thermal)
    E⁻[:,1]=boundary_top;
    # Source at bottom layer (0 for SIF, emitted from soil for thermal)
    U[:,end].=boundary_bottom;

    @inbounds for j=nl:-1:1   # bottom to top
        dnorm       = 1 .-ρ_dd[:,j].*Rdd[:,j+1];
        Xdd[:,j]  = τ_dd[:,j]./dnorm;
        Rdd[:,j]  = ρ_dd[:,j]+τ_dd[:,j].*Rdd[:,j+1].*Xdd[:,j];
        Y[:,j] = (ρ_dd[:,j].*U[:,j+1]+S⁻[:,j])./dnorm;
        U[:,j] = τ_dd[:,j].*(Rdd[:,j+1].*Y[:,j]+U[:,j+1])+S⁺[:,j];
    end
    @inbounds for j=1:nl      # from top to bottom
        E⁻[:,j+1] = Xdd[:,j].*E⁻[:,j]+Y[:,j];
        E⁺[:,j]   = Rdd[:,j].*E⁻[:,j]+U[:,j];
    end
    E⁺[:,end] = Rdd[:,end].*E⁻[:,end]+U[:,end];
    # compute net diffuse radiation per layer:
    @inbounds for j=1:nl
        net_diffuse[:,j]= (E⁻[:,j] +  E⁺[:,j+1]).*(1 .-(τ_dd[:,j]+ρ_dd[:,j]))
    end
    return E⁻,E⁺,net_diffuse

end;


# Postprocessing, computing all within canopy fluxes from the primary RT output
function deriveCanopyFluxes!(can::struct_canopy, cO::struct_canopyOptProps, cR::struct_canopyRadiation, sun::incomingRadiation, so::struct_soil, leaf::Array)
    #
    nl = can.nlayers
    iLAI    = (can.LAI*can.Ω)/nl;
    rsoil = so.albedo_SW;
    soil_emissivity = 1 .-rsoil # emissivity of soil (goes into structure later)


    # Compute some fluxes, can be done separately if needed (this is absolute fluxes now, for the entire soil):
    cR.RnSoil_diffuse = fac * fast∫(dwl, cR.E_down[:,end].*soil_emissivity);
    cR.RnSoil_direct  = fac *fast∫(dwl, cO.Es_[:,end].*soil_emissivity);
    cR.RnSoil = cR.RnSoil_direct + cR.RnSoil_diffuse

    # Normalization factor for leaf direct PAR (weighted sum has to be 1 to conserve net SW direct)
    normi = 1/mean(cO.absfs'*can.lidf)
    lPs = (cO.Ps[1:nl]+cO.Ps[2:nl+1])/2
    @inbounds for j = 1:nl
        if length(leaf)>1
            	kChlrel = leaf[j].kChlrel[iPAR]
        else
        kChlrel = leaf[1].kChlrel[iPAR]
        end
        PAR_diff = (fac/iLAI)*e2phot(wl[iPAR], cR.netSW_shade[iPAR,j])
        # Direct PAR is normalized by layer Ps value:
        PAR_dir = (fac/iLAI/lPs[j])*e2phot(wl[iPAR], cR.netSW_sunlit[iPAR,j])+PAR_diff

        PAR_diffCab = kChlrel.*PAR_diff
        #println(leaf.kChlrel[iPAR])
        # Direct PAR is normalized by layer Ps value:
        PAR_dirCab  = kChlrel.*PAR_dir;

        # Absorbed PAR per leaf for shaded:
        cR.absPAR_shade[j] = fast∫(dwl[iPAR],PAR_diff);
        # for sunlit (all angles, normalized)
        #@show normi*cO.absfs*fast∫(dwl[iPAR],PAR_dir)
        cR.absPAR_sun[:,:,j]  = normi*cO.absfs*fast∫(dwl[iPAR],PAR_dir);
        # Same for true absorbed PAR by CabCar only  (needs to be taken into account in photosynthesis calculation, i.e. already accounts for leaf green fAPAR!)
        cR.absPAR_shadeCab[j] = fast∫(dwl[iPAR],PAR_diffCab);
        cR.absPAR_sunCab[:,:,j]  = normi*cO.absfs*fast∫(dwl[iPAR],PAR_dirCab);
    #println(PAR_dir[1], " ", cR.Pnh[j])
    end
    #@time fast∫(dwl[iPAR], sun.E_direct[iPAR])
    cR.incomingPAR_direct   = fac * fast∫(dwl[iPAR], e2phot(wl[iPAR],(sun.E_direct[iPAR])));
    cR.incomingPAR_diffuse  = fac * fast∫(dwl[iPAR], e2phot(wl[iPAR],(sun.E_diffuse[iPAR])));
    cR.incomingPAR          = cR.incomingPAR_diffuse*cR.incomingPAR_direct
    @inbounds for i=1:nl
        cR.intNetSW_shade[i] =   (fac/iLAI) *fast∫(dwl, cR.netSW_shade[:,i]);
        cR.intNetSW_sunlit[i]  =   (fac/iLAI/lPs[i]) * fast∫(dwl, cR.netSW_sunlit[:,i]) + cR.intNetSW_shade[i];
    end

end;

function compCanopyOptsExact!(leaf::Array,cO::struct_canopyOptProps, lidf::Array)
    @unpack fs,fo, fsfo, absfsfo, cosΘ_l = cO
    fill!(cO.w,0)
    fill!(cO.sf,0)
    fill!(cO.sb,0)
    fill!(cO.vf,0)
    fill!(cO.vb,0)
    fill!(cO.sigf,0)
    fill!(cO.sigb,0)


    f₁ = lidf.*(1 .+ cosd.(litab))/2
    f₂ =  lidf.*(1 .- cosd.(litab))/2
    @inbounds for i=1:nlayers
        if length(leaf)>1
            τ = leaf[i].τ_SW
            ρ = leaf[i].ρ_SW
        else
            τ = leaf[1].τ_SW
            ρ = leaf[1].ρ_SW
        end
        R = CartesianIndices(canOpt.fs)
        for I in R
            iL = I[1]
            @inbounds for j=1:size(cO.sigb, 1)
                    cO.w[j,i]  += lidf[iL]*(abs(fsfo[I])*(ρ[j]+τ[j])+fsfo[I]*(ρ[j]-τ[j]))/2
                    cO.vb[j,i] += lidf[iL]*(abs(fo[I])*(ρ[j]+τ[j])+fo[I]*cosΘ_l[I]*(ρ[j]-τ[j]))/2; # [nl,nwl]  directional backscatter scattering coefficient for diffuse  incidence
                    cO.vf[j,i] += lidf[iL]*(abs(fo[I])*(ρ[j]+τ[j])-fo[I]*cosΘ_l[I]*(ρ[j]-τ[j]))/2 # [nl,nwl]  directional forward     scattering coefficient for diffuse  incidence
                    #cO.w[:,i] += lidf[iL]*(abs(fsfo[I])*(ρ.+τ)+fsfo[I]*(ρ.-τ))/2
                    #w2[:] += canopy.lidf[I[1]]*(abs(canOpt.fsfo[I])*(rho.+tau)+canOpt.fsfo[I]*(rho.-tau))/2
            end
        end
    end
end;


function computeCanopyMatrices!(leaf::Array,cO::struct_canopyOptProps)
    # 2. Calculation of reflectance
    # 2.1  reflectance, transmittance factors in a thin layer the following are vectors with length [nl,nwl]
    @inbounds for i=1:nlayers
    if length(leaf)>1
        	τ_SW = leaf[i].τ_SW
        	ρ_SW = leaf[i].ρ_SW
    else
    τ_SW = leaf[1].τ_SW
        	ρ_SW = leaf[1].ρ_SW
    end
    #CF: Right now, canopy geometry is the same everywhere, can be easily extended to layers as well.
        @inbounds for j=1:size(cO.sigb, 1)
          cO.sigb[j,i]= cO.ddb  * ρ_SW[j]   + cO.ddf * τ_SW[j]; # [nl,nwl]  diffuse     backscatter scattering coefficient for diffuse  incidence
          cO.sigf[j,i]= cO.ddf  * ρ_SW[j]   + cO.ddb * τ_SW[j]; # [nl,nwl]  diffuse     forward     scattering coefficient for diffuse  incidence
          cO.sb[j,i]  = cO.sdb  * ρ_SW[j]   + cO.sdf * τ_SW[j]; # [nl,nwl]  diffuse     backscatter scattering coefficient for specular incidence
          cO.sf[j,i]  = cO.sdf  * ρ_SW[j]   + cO.sdb * τ_SW[j]; # [nl,nwl]  diffuse     forward     scattering coefficient for specular incidence
          cO.vb[j,i]  = cO.dob  * ρ_SW[j]   + cO.dof * τ_SW[j]; # [nl,nwl]  directional backscatter scattering coefficient for diffuse  incidence
          cO.vf[j,i]  = cO.dof  * ρ_SW[j]   + cO.dob * τ_SW[j]; # [nl,nwl]  directional forward     scattering coefficient for diffuse  incidence
          cO.w[j,i]   = cO.sob  * ρ_SW[j]   + cO.sof * τ_SW[j]; # [nl,nwl]  bidirectional scattering coefficent (directional-directional)

        end
    end
    cO.a    .= 1 .- cO.sigf;                  # [nl, nwl]     attenuation
end;

# Compute optical properties of canopy
function computeCanopyGeomProps!(can::struct_canopy, angle::struct_angles, cO::struct_canopyOptProps)
    @unpack Ω,nlayers,LAI = can
    # Number of layers
    nl = nlayers

    dx = FT(1/nl)
    xl = collect(FT(0):FT(-1/nl):-1)
    tts = angle.tts
    tto = angle.tto
    psi = angle.psi
    lazitab = can.lazitab
    # only needed for volume scattering for symmetry (not sure why it wasn't working)
    psi_vol         = abs(psi-FT(360.0)*round(psi/FT(360.0)));
    # Geometric quantities (ougoing direction first!)
    cto = cosd(tto)
    tanto	= tand(tto);
    cospsi	= cosd(psi);

    # Generate leaf angle distribution:
    
    lidf = dladgen(can.LIDFa,can.LIDFb, litab_bnd);
    
    can.lidf = lidf

    #println(lidf)

    # Precompute all of this before in a separate function?
    cos_ttlo    = cosd.(lazitab);  #   cos leaf azimuth angles
    cos_philo    = cosd.(lazitab.-psi);  #   cos leaf azimuth angles
    cos_ttli    = cosd.(litab);    #   cosine of normal of upperside of leaf
    sin_ttli    = sind.(litab);    #   sine   of normal of upperside of leaf
    sin_tto  = sind(tto);
    # Solar angle dependent ones:
    cts = cosd(tts)
    sin_tts  = sind(tts);
    ctscto = cts*cto;
    tants	= tand(tts);
    dso		= sqrt(tants*tants+tanto*tanto-2tants*tanto*cospsi);

    # angular distance, compensation of shadow length
    # Calculate geometric factors associated with extinction and scattering
    #Initialise sums
    cO.ks    = 0;
    cO.ko    = 0;
    cO.bf    = 0;
    cO.sob = 0;
    cO.sof  = 0;

    #	Weighted sums over LIDF
    @inbounds for i=1:length(litab)
        # ttl = litab[i];	% leaf inclination discrete values
        ctl = cosd(litab[i]);
        #	SAIL volume scattering phase function gives interception and portions to be multiplied by rho and tau
        chi_s,chi_o,frho,ftau=volscatt(tts,tto,psi_vol,litab[i]);
        #	Extinction coefficients
        ksli = abs(chi_s./cts);
        koli =abs(chi_o./cto);
            #	Area scattering coefficient fractions
        #@show ftau
        #@show frho
        if ftau<0
            sobli	= abs(ftau)*FT(pi)/ctscto;
            sofli	= abs(frho)*FT(pi)/ctscto;
        else
            sobli	=(frho*FT(pi)/ctscto);
            sofli	= (ftau*FT(pi)/ctscto);
        end
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
    cds  = Cs*ones(FT,1,length(lazitab)) .+ Ss*cos_ttlo';                    # [nli,nlazi]
    cdo  = (cos_ttli.*cto)*ones(FT,1,length(lazitab)) .+ (sin_ttli.*sin_tto)*cos_philo';   # [nli,nlazi]

    # This is basically equivalent to Kb in Bonan, eq. 14.21
    cO.fs          = (cds./cts);
    cO.absfs    = abs.(cO.fs);         # [nli,nlazi] pag 305
    # Added here
    cO.fo = cdo/cto;
    cO.cosΘ_l = cos_ttli.*ones(FT,1,length(lazitab));
    cO.cos2Θ_l = cO.cosΘ_l.^2;
    cO.fsfo = (cO.fs.*cO.fo)
    cO.absfsfo = abs.(cO.fsfo)

    # 1.5 probabilities Ps, Po, Pso
    cO.Ps[:]          =   exp.(cO.ks*xl*Ω*LAI);     # [nl+1]  probability of viewing a leaf in solar dir
    # Leave off Omega for Po here still
    cO.Po[:]          =   exp.(cO.ko*xl*Ω*LAI);     # [nl+1]  probability of viewing a leaf in observation dir

    cO.Ps    =   cO.Ps *(1 .-exp.(-cO.ks*Ω*LAI*dx))/(cO.ks*Ω*LAI*dx);   # Correct Ps/Po for finite dx
    cO.Po    =   cO.Po *(1 .-exp.(-cO.ko*Ω*LAI*dx))/(cO.ko*Ω*LAI*dx);   # Correct Ps/Po for finite dx

    #Pso: Probability of observing a sunlit leaf at depth x, see eq 31 in vdT 2009
    # Just ignore Ω for now, it is not yet thought through!
    #@show cO.ks
    #@show cO.ko
    @inbounds for j=1:length(xl)
        cO.Pso[j] = quadgk(x -> Psofunction(cO.ko,cO.ks,Ω,LAI,can.hot,dso,x), xl[j]-dx,xl[j], rtol=1e-2)[1]/dx
    end

    #cO.Pso[cO.Pso.>cO.Po]= minimum([cO.Po[cO.Pso.>cO.Po] cO.Ps[cO.Pso.>cO.Po]],dims=2);    #takes care of rounding error
    #cO.Pso[cO.Pso.>cO.Ps]= minimum([cO.Po[cO.Pso.>cO.Ps] cO.Ps[cO.Pso.>cO.Ps]],dims=2);    #takes care of rounding error
end;

function computeThermalFluxes!( leaf::Array,cO::struct_canopyOptProps,cR::struct_canopyRadiation, can::struct_canopy,so::struct_soil, incLW::Array)

    @unpack Ps, Po, Pso,ddf,ddb = cO
    @unpack T_sun, T_shade = cR
    @unpack Ω,nlayers, LAI = can
    @unpack albedo_LW, soil_skinT = so
    FT = eltype(LAI)

    iLAI    = Ω*LAI/nlayers;

    # Sunlit fraction at mid layer:
    fSun      = (Ps[1:nlayers] + Ps[2:nlayers+1])/2;

    ϵ = similar(T_shade)
    τ_dd = zeros(FT,1,length(T_shade))
    ρ_dd = similar(τ_dd)
    S⁺ = similar(ρ_dd)
    S⁻ = similar(ρ_dd)
    S_shade = similar(ρ_dd)
    S_sun = similar(ρ_dd)


    # If we just have the same leaf everywhere, compute emissivities:
    if length(leaf)==1
        le = leaf[1]
        # Compute layer properties:
        sigf = ddf*le.ρ_LW + ddb*le.τ_LW
        sigb = ddb*le.ρ_LW + ddf*le.τ_LW
        τ_dd = (1 - (1-sigf)*iLAI)*ones(nwl,nl)
        ρ_dd = (sigb*iLAI)*ones(nwl,nl)
        ϵ .= (1 - τ_dd-ρ_dd);
    elseif length(leaf)==nlayers
        for i=1:nlayers
            le = leaf[i]
            sigf = ddf*le.ρ_LW + ddb*le.τ_LW
            sigb = ddb*le.ρ_LW + ddf*le.τ_LW
            τ_dd[i] = (1 - (1-sigf)*iLAI)
            #τ_dd[i] = (1 - exp(-sigf*iLAI))
            ρ_dd[i]= (sigb*iLAI)
            ϵ[i] = (1 - τ_dd[i]-ρ_dd[i]);
        end
    else
        println("Complain, Array of leaves is neither 1 nor nlayers ")
    end

    # Only one wavelength --> do Stefan Boltzmann:
    #if length(wl)==1
    # Let's just do SB for now:
    if 1==1
        # Shaded leaves first, simple 1D array:
        S_shade=physcon.σ.*ϵ.*(T_shade.^4)
        # Sunlit leaves:
        if ndims(T_sun)>1
            @inbounds for i=1:length(T_shade)
                emi = physcon.σ*ϵ[i]*T_sun[:,:,i].^4
                # weighted average over angular distribution
                S_sun[i] = mean(emi'*lidf);
            end
        else
            # Sunlit, simple 1D array:
            S_sun = physcon.σ.*ϵ.*T_sun.^4
        end
    else
        # Do Planck curve, tbd
    end
    S⁺[:] = iLAI*(fSun.*S_sun+(1 .-fSun).*S_shade)
    S⁻[:] = S⁺[:]
    soilEmission = physcon.σ*(1 .- albedo_LW) * soil_skinT^4
    # Run RT:
    F⁻,F⁺,net_diffuse = RTM_diffuseS(τ_dd, ρ_dd,S⁻, S⁺,incLW, soilEmission, albedo_LW)
    for j = 1:nlayers
        cR.intNetLW_sunlit[j]   = net_diffuse[j]- 2S_sun[j];          #  sunlit leaf
        cR.intNetLW_shade[j]   = net_diffuse[j]- 2S_shade[j];     # shaded leaf
    end
    # Net soil LW as difference between up and downwelling at lowest level
    cR.RnSoilLW = F⁻[1,end]-F⁺[1,end]
    #@show F⁻[:,end]
    #@show F⁺[:,end]
    return F⁻,F⁺,net_diffuse
end;

function computeSIF_Fluxes!(leaf::Array,cO::struct_canopyOptProps,cR::struct_canopyRadiation,can::struct_canopy,so::struct_soil)
    @unpack fs,fo, cosΘ_l, absfs, fsfo, absfsfo, cos2Θ_l, Ps, Po, Pso, a,  sigb,vb,vf = cO
    @unpack E_down, E_up,φ_shade,φ_sun = cR
    @unpack Ω,nlayers, LAI = can
    rsoil = so.albedo_SW[Iwlf]
    iLAI    = Ω*LAI/nlayers;

    τ_dd = 1 .-a[Iwlf,:]*iLAI;
    ρ_dd = sigb[Iwlf,:]*iLAI;

    dim = size(leaf[1].Mb)
    # Compute mid layer Ps,Po,Pso
    #Qso   = (Pso[1:nlayers] + Pso[2:nlayers+1])/2;
    #Qs      = (Ps[1:nlayers] + Ps[2:nlayers+1])/2;
    #Qo      =(Po[1:nlayers] + Po[2:nlayers+1])/2;
    Qso   = Pso[1:nlayers]
    Qs     = Ps[1:nlayers]
    Qo     = Po[1:nlayers]
    # Allocate arrays for output:
    di = size(leaf[1].Mb)
    S⁻ = zeros(FT,dim[1],nlayers )
    S⁺  = similar(S⁻)
    piLs = similar(S⁻)
    piLd = similar(S⁻)
    Fsmin = similar(S⁻)
    Fsplu = similar(S⁻)
    Fdmin = similar(S⁻)
    Fdplu = similar(S⁻)
    Femo = similar(S⁻);

    # fs is different here than in RTMo (in RTMo, it is abs(fs))!

    # 2. Calculation of reflectance
    # 2.1  reflectance, transmittance factors in a thin layer the following are vectors with length [nl,nwl]
    sun = cO.Es_[Iwle,1];
    @inbounds for i=1:nlayers
        if length(leaf)>1
            Mb = leaf[i].Mb
            Mf = leaf[i].Mf
        else
            Mb = leaf[1].Mb
            Mf = leaf[1].Mf
        end
        M⁺ = (Mb + Mf)/2;
        M⁻ = (Mb - Mf)/2;
        # Need to normalize incoming radiation bin, change from mSCOPE to enable different wl grids!
        M⁻_sun = M⁻ * (sun.*dwl[Iwle]);
        M⁺_sun = M⁺ * (sun.*dwl[Iwle]);

        M⁺⁻ = M⁺ * (E_down[Iwle,i].*dwl[Iwle]);
        M⁺⁺ = M⁺ * (E_up[Iwle,i+1].*dwl[Iwle]);
        M⁻⁺ = M⁻ * (E_up[Iwle,i+1].*dwl[Iwle]);
        M⁻⁻ = M⁻ * (E_down[Iwle,i].*dwl[Iwle]);


        # Here comes the tedious part:
        sunCos    = mean((φ_sun[:,:,i].*cosΘ_l)'*can.lidf)
        shadeCos  = mean((φ_shade[i]*cosΘ_l)'*can.lidf)
        sunCos2   = mean((φ_sun[:,:,i].*cos2Θ_l)'*can.lidf)
        shadeCos2 = mean((φ_shade[i]*cos2Θ_l)'*can.lidf)
        sunLidf   = mean(φ_sun[:,:,i]'*can.lidf)
        shadeLidf = mean(φ_shade[i]'*can.lidf)

        wfEs = mean((φ_sun[:,:,i].*absfsfo)'*can.lidf)  * M⁺_sun + mean((φ_sun[:,:,i].*fsfo)'*can.lidf)  *M⁻_sun

        a1 = mean((φ_sun[:,:,i].*absfs)'*can.lidf)* M⁺_sun;
        a2 = mean((φ_sun[:,:,i].*fs.*cosΘ_l)'*can.lidf)*M⁻_sun
        sfEs = a1 - a2
        sbEs = a1 + a2

        a1 = mean((φ_shade[i]*abs.(fo))'*can.lidf) ;
        a2 = mean((φ_shade[i]*fo.*cosΘ_l)'*can.lidf)
        vfEplu_shade =  a1  * M⁺⁺ + a2 *M⁻⁺
        vbEmin_shade =  a1 * M⁺⁻ + a2 *M⁻⁻

        a1 = mean((φ_sun[:,:,i].*abs.(fo))'*can.lidf);
        a2 = mean((φ_sun[:,:,i].*fo.*cosΘ_l)'*can.lidf)
        vfEplu_sun  =  a1 * M⁺⁺ - a2 *M⁻⁺
        vbEmin_sun  =  a1 * M⁺⁻ + a2 *M⁻⁻

        a1 = shadeLidf * M⁺⁻;
        a2 = shadeCos2 * M⁻⁻
        sigfEmin_shade  =  a1 - a2
        sigbEmin_shade  =  a1 + a2

        a1 = sunLidf * M⁺⁻;
        a2 = sunCos2 * M⁻⁻;
        sigfEmin_sun  =  a1 - a2
        sigbEmin_sun  =  a1 +  a2

        a1 = shadeLidf  * M⁺⁺;
        a2 = shadeCos2  * M⁻⁺
        sigfEplu_shade  = a1- a2
        sigbEplu_shade  = a1 + a2

        a1 = sunLidf  * M⁺⁺;
        a2 = sunCos2  * M⁻⁺;
        sigfEplu_sun  = a1- a2
        sigbEplu_sun  = a1 +a2

        # Fluxes:
        piLs[ :,i]  =  wfEs+vfEplu_sun+vbEmin_sun;       # sunlit for each layer
        piLd[ :,i]  =  vbEmin_shade+vfEplu_shade;        # shade leaf for each layer
        Fsmin[:,i]  =  sfEs+sigfEmin_sun+sigbEplu_sun;   # Eq. 29a for sunlit leaf
        Fsplu[:,i]  =  sbEs+sigbEmin_sun+sigfEplu_sun;   # Eq. 29b for sunlit leaf
        Fdmin[:,i]  =  sigfEmin_shade+sigbEplu_shade;    # Eq. 29a for shade leaf
        Fdplu[:,i]  =  sigbEmin_shade+sigfEplu_shade;    # Eq. 29b for shade leaf
        # Total weighted fluxes
        S⁻[:,i]     =   iLAI*(Qs[i]*Fsmin[:,i] + (1-Qs[i])*Fdmin[:,i]);
        S⁺[:,i]     =   iLAI*(Qs[i]*Fsplu[:,i] + (1-Qs[i])*Fdplu[:,i]);
        
        Femo[:,i]   =   iLAI*(Qs[i]* piLs[:,i] + (1-Qs[i])*piLd[:,i]);
    end
    # Use Zero SIF fluxes as top and bottom boundary:
    zeroB = zeros(FT,length(Iwlf))
    # Compute diffusive fluxes within canopy
    #println(size(zeroB), " ", size(rsoil), " ", size(S⁻), " ", size(τ_dd))
    F⁻,F⁺,net_diffuse = RTM_diffuseS(τ_dd, ρ_dd,S⁻, S⁺,zeroB, zeroB, rsoil)
    # Save in output structures!
    cR.SIF_obs_sunlit[:]  = iLAI/FT(pi)*Qso'*piLs';                                               # direct Sunlit leaves
    cR.SIF_obs_shaded[:]  = iLAI/FT(pi)*(Qo[1:nlayers]-Qso[1:nlayers])'*piLd';                    # direct shaded leaves

    # SIF scattered internally
    #cR.SIF_obs_scattered[:]     = iLAI/FT(pi)*(Qo[1:nlayers]'*(vb[Iwlf,:].*F⁻[:,1:nlayers] + vf[Iwlf,:].*F⁺[:,1:nlayers])');
    cR.SIF_obs_scattered[:]     = iLAI/FT(pi)*(Qo[1:nlayers]'*(vb[Iwlf,:].*F⁻[:,1:nlayers] + vf[Iwlf,:].*F⁺[:,1:nlayers])');
    cR.SIF_obs_soil[:]     = (rsoil .* F⁻[:,end] * Po[end])/FT(pi);                                                #Soil contribution
    
    cR.SIF_hemi[:]    = F⁺[:,1];
    cR.SIF_obs[:] = cR.SIF_obs_sunlit[:]+cR.SIF_obs_shaded[:]+cR.SIF_obs_scattered[:] +cR.SIF_obs_soil[:];
    cR.SIF_sum[:] = sum(S⁻+S⁺, dims=2)
    return
    #return  F⁻,F⁺,S⁻,S⁺, piLs, piLd
end;


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
    @unpack N, Cab, Car, Ant, Cs, Cw, Cm, ρ_SW, τ_SW, Cx = leaf
    #println(N, " ", Cab, " ", Car," ",  Ant, " ", Cs, " ", Cw, " ", Cm)
    Kcaro = (1 -Cx)* opti.KcaV + Cx * opti.KcaZ;

    Kall    = (Cab*opti.Kab.+Car*opti.Kcar.+Ant*opti.Kant.+Cs*opti.KBrown.+Cw*opti.Kw.+Cm*opti.Km)/N
    # Relative absorption by Chlorophyll and Carotenoids only (drives SIF and GPP eventually)
    leaf.kChlrel  = (Cab*opti.Kab+Car*opti.Kcar)./(Kall*N.+eps(FT));
    leaf.kChlrel_old  = (Cab*opti.Kab)./(Kall*N.+eps(FT));
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
    ralf    = FT(1) .-talf

    t12     = calctav.(90,opti.nr)
    r12     = FT(1) .-t12
    t21     = t12./(opti.nr.^2)
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
    Rsub[j]	= 1 .-Tsub[j]

    # Reflectance & transmittance of the leaf: combine top layer with next N-1 layers
    denom   = 1 .-Rsub.*r
    leaf.τ_SW    = Ta.*Tsub./denom
    leaf.ρ_SW    = Ra.+Ta.*Rsub.*t./denom
    τ_SW    = leaf.τ_SW
    ρ_SW    = leaf.ρ_SW
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

    # indices of wle and wlf within wlp

    epsi         = FT(2)^(-ndub);

    # initialisations
    te     = 1 .-(k[Iwle].+s[Iwle]) * epsi;
    tf     = 1 .-(k[Iwlf].+s[Iwlf]) * epsi;
    re     = s[Iwle] * epsi;
    rf     = s[Iwlf] * epsi;

    sigmoid     = 1 ./(1 .+exp.(-wlf/10).*exp.(wle'/10));  # matrix computed as an outproduct
    #println(size(sigmoid)," ", size(phi), " ", size(kChl)," ", size(Iwle), " ", size(Iwlf), " ", size(kChl[Iwle]))

    Mf = leaf.fqe .* ((FT(0.5)*opti.phi[Iwlf]).*epsi) .* kChl[Iwle]'.*sigmoid
    Mb = leaf.fqe .* ((FT(0.5)*opti.phi[Iwlf]).*epsi) .* kChl[Iwle]'.*sigmoid

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
    # Here we add the leaf-air interfaces again for obtaining the final
    # leaf level fluorescences.

    # This reduced red SIF quite a bit in backscatter, not sure why.
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
    return Mf,Mb
end;



# Had to make this like the Prosail F90 function, forwardDiff didn't work otherwise.
function calcJ1(x,m,k,LAI,Ω)
    del=(k-m)*LAI
    if abs(del)>1E-3;
        J1 = (exp(m*Ω*LAI*x)-exp(k*Ω*LAI*x))/(k-m);
    else
        J1 = -0.5*Ω*LAI*x*(exp(m*Ω*LAI*x)+exp(k*Ω*LAI*x))*(1.0-del^2/12.0);
    end
    return J1
end

function calcJ2(x,m,k,LAI,Ω)
    return (exp(k*Ω*LAI*x)-exp(-k*Ω*LAI)*exp(-m*Ω*LAI*(1+x)))/(k+m);
end




# APPENDIX IV function Pso from SCOPE v1.73
function Psofunction(K,k,Ω,LAI,q,dso,xl)
    if dso!=0.0
        alf         =   (dso/q) *2/(k+K);
        pso         =   exp((K+k)*Ω*LAI*xl + sqrt(K*k)*Ω*LAI/(alf  )*(1-exp(xl*(alf  ))));# [nl+1]  factor for correlation of Ps and Po
    else
        pso         =   exp((K+k)*Ω*LAI*xl - sqrt(K*k)*Ω*LAI*xl);# [nl+1]  factor for correlation of Ps and Po
    end
    return pso
end;

# Changed from latest volscat functions (April 2001)
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
    #@show Ss
    #@show So
    #@show Cs
    #@show Co
    cosbts=1;
    if (abs(Ss)>1e-6)
    	cosbts=-Cs/Ss;
    end
    cosbto=1;
    if (abs(So)>1e-6)
    	cosbto=-Co/So;
    end

    if (abs(cosbts)<1)
    	bts=acos(cosbts);
    	ds=Ss;
    else
    	bts=FT(pi);
    	ds=Cs;
    end

    chi_s=2/FT(pi)*((bts-FT(pi)/2)*Cs+sin(bts)*Ss);

    if (abs(cosbto)<1)
    	bto=acos(cosbto);
    	doo=So;
    elseif(tto<90)
    	bto=pi;
    	doo=Co;
    else
    	bto=0;
    	doo=-Co;
    end
    chi_o=2 /FT(pi)*((bto-FT(pi)/2)*Co+sin(bto)*So);

    #c ..............................................................................
    #c   Computation of auxiliary azimut angles bt1, bt2, bt3 used
    #c   for the computation of the bidirectional scattering coefficient w
    #c .............................................................................

    btran1=abs(bts-bto);
    btran_=pi-abs(bts+bto-FT(pi));
    btran2=2FT(pi)-bts-bto;
    #@show psi_rad
    #@show btran1
    #@show btran2

    if (psi_rad<=btran1)
    	bt1=psi_rad;
    	bt2=btran1;
    	bt3=btran2;
    else
    	bt1=btran1;
    	if (psi_rad<=btran2)
    		bt2=psi_rad;
    		bt3=btran2;
        else
    		bt2=btran2;
    		bt3=psi_rad;
        end
    end

    t1=2Cs*Co+Ss*So*cos_psi;
    t2=0;
    if (bt2>0)
    	t2=sin(bt2)*(2ds*doo+Ss*So*cos(bt1)*cos(bt3));
    end

    denom=2 *FT(pi)^2;
    frho=((pi-bt2)*t1+t2)/denom;
    ftau=    (-bt2*t1+t2)/denom;
    #if frho<0 || ftau<0 || chi_s<0 || chi_o < 0
        #@show bto
        #@show bts
        #@show bts-bto
        #@show btran1
        #@show btran2
        #@show btran_
        #@show psi_rad
    #end
    return (chi_s),abs(chi_o),(frho),(ftau)
end;

function dladgen(a::Number,b::Number, litab::Array)
    @assert a+b<1
    freq = similar(litab[:,1])
    for i in eachindex(freq)
        freq[i]=dcum(a,b,litab[i,2])-dcum(a,b,litab[i,1]);
    end
    return freq
end;

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
            y  = a*sin(x)+0.5*b*sin(2x);
            dx = (y-x+p)/2;
            x=x+dx;
            delx=abs(dx);
        end
    	f = (2.0*y+p)/pi;
    end
    return f
end;


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
