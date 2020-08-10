###############################################################################
#
# Update canopy geometry
#
###############################################################################
"""
    canopy_geometry!(can::Canopy4RT{FT}, angles::SolarAngles{FT}, can_opt::CanopyOpticals{FT}) where {FT<:AbstractFloat}

Computes canopy optical properties (extinction coefficients for direct and
    diffuse light) based on the SAIL model. Most important input parameters are
    leaf inclination and azimuth distribution functions and sun-sensor
    geometry . Canopy clumping Ω is implemented as in Pinty et al (2015), given
- `can` A [`Canopy4RT`](@ref) type struct, providing canopy structure information
- `angles` A [`SolarAngles`](@ref) type struct, defining sun-sensor geometry
- `can_opt` An [`CanopyOpticals`](@ref) type, where optical properties will be stored

"""
function canopy_geometry!(
            can::Canopy4RT{FT},
            angles::SolarAngles{FT},
            can_opt::CanopyOpticals{FT}
) where {FT<:AbstractFloat}
    @unpack clump_a, clump_b, LAI, litab, litab_bnd, nLayer, Ω,lazitab = can
    @unpack tts,tto,psi = angles

    if clump_b > 0
        can.Ω = clump_a .+ clump_b.*(1. .- cosd(tts))
    end

    dx  = FT(1/nLayer)
    xl  = collect(FT(0):FT(-1/nLayer):-1)

    # only needed for volume scattering for symmetry (not sure why it wasn't working)
    psi_vol = abs(psi-FT(360.0)*round(psi/FT(360.0)))
    # Geometric quantities (ougoing direction first!)
    cto     = cosd(tto)
    tanto   = tand(tto)
    cospsi  = cosd(psi)

    # Generate leaf angle distribution:
    # This is one-time calculation, move it to when initializing the struct
    lidf     = dladgen(can.LIDFa, can.LIDFb, litab_bnd)
    can.lidf = lidf

    # Precompute all of this before in a separate function?
    cos_ttlo   = cosd.(lazitab)           # cos leaf azimuth angles
    cos_philo  = cosd.(lazitab .- psi)    # cos leaf azimuth angles
    cos_ttli   = cosd.(litab)             # cos of normal of upperside of leaf
    sin_ttli   = sind.(litab)             # sine of normal of upperside of leaf
    sin_tto    = sind(tto)
    # Solar angle dependent ones:
    cts        = cosd(tts)
    sin_tts    = sind(tts)
    ctscto     = cts*cto
    tants	   = tand(tts)
    dso		   = sqrt(tants*tants+tanto*tanto-2tants*tanto*cospsi)

    # Calculate geometric factors associated with extinction and scattering
    can_opt.ks  = 0
    can_opt.ko  = 0
    can_opt.bf  = 0
    can_opt.sob = 0
    can_opt.sof = 0

    #	Weighted sums over LIDF
    @inbounds for i=1:length(litab)
        # ttl = litab[i]	% leaf inclination discrete values
        ctl = cosd(litab[i])
        # SAIL volume scattering phase function gives interception and portions
        # to be multiplied by rho and tau
        chi_s,chi_o,frho,ftau=volscatt(tts,tto,psi_vol,litab[i])
        # Extinction coefficients
        ksli = abs(chi_s./cts)
        koli = abs(chi_o./cto)
        # Area scattering coefficient fractions
        #@show ftau
        #@show frho
        if ftau<0
            sobli = abs(ftau)*FT(pi)/ctscto
            sofli = abs(frho)*FT(pi)/ctscto
        else
            sobli = (frho*FT(pi)/ctscto)
            sofli = (ftau*FT(pi)/ctscto)
        end
        bfli   = ctl*ctl
        can_opt.ks  = can_opt.ks+ksli*lidf[i]
        can_opt.ko  = can_opt.ko+koli*lidf[i]
        can_opt.bf  = can_opt.bf+bfli*lidf[i]
        can_opt.sob = can_opt.sob+sobli*lidf[i]
        can_opt.sof = can_opt.sof+sofli*lidf[i]
    end
    #println(sob, " ", sof)

    #	Geometric factors to be used later with rho and tau
    can_opt.sdb = (can_opt.ks+can_opt.bf)/2
    can_opt.sdf = (can_opt.ks-can_opt.bf)/2
    can_opt.dob = (can_opt.ko+can_opt.bf)/2
    can_opt.dof = (can_opt.ko-can_opt.bf)/2
    can_opt.ddb = (1 .+can_opt.bf)/2
    can_opt.ddf = (1 .-can_opt.bf)/2

    # See eq 19 in vdT 2009
    Cs  = cos_ttli.*cts             # [nli]     pag 305 modified by Joris
    Ss  = sin_ttli.*sin_tts         # [nli]     pag 305 modified by Joris
    cds = Cs*ones(FT,1,length(lazitab)) .+ Ss*cos_ttlo'   # [nli,nlazi]
    cdo = (cos_ttli.*cto)*ones(FT,1,length(lazitab)) .+
          (sin_ttli.*sin_tto)*cos_philo'                  # [nli,nlazi]

    # This is basically equivalent to Kb in Bonan, eq. 14.21
    can_opt.fs      = (cds./cts)
    can_opt.absfs   = abs.(can_opt.fs)         # [nli,nlazi] pag 305
    # Added here
    can_opt.fo      = cdo/cto
    can_opt.cosΘ_l  = cos_ttli.*ones(FT,1,length(lazitab))
    can_opt.cos2Θ_l = can_opt.cosΘ_l.^2
    can_opt.fsfo    = (can_opt.fs.*can_opt.fo)
    can_opt.absfsfo = abs.(can_opt.fsfo)

    # 1.5 probabilities Ps, Po, Pso
    # probability of viewing a leaf in solar dir
    can_opt.Ps[:] = exp.(can_opt.ks*xl*Ω*LAI)     # [nl+1]
    # Leave off Omega for Po here still
    # probability of viewing a leaf in observation dir
    can_opt.Po[:] = exp.(can_opt.ko*xl*Ω*LAI)     # [nl+1]

    # # Correct Ps/Po for finite dx
    can_opt.Ps = can_opt.Ps *(1 .-exp.(-can_opt.ks*Ω*LAI*dx))/(can_opt.ks*Ω*LAI*dx)
    can_opt.Po = can_opt.Po *(1 .-exp.(-can_opt.ko*Ω*LAI*dx))/(can_opt.ko*Ω*LAI*dx)

    #Pso: Probability of observing a sunlit leaf at depth x, see eq 31 in vdT 2009
    # Just ignore Ω for now, it is not yet thought through!
    #@show can_opt.ks
    #@show can_opt.ko
    @inbounds for j=1:length(xl)
        can_opt.Pso[j] = quadgk(x -> psofunction(can_opt.ko,can_opt.ks,Ω,LAI,can.hot,dso,x),
                           xl[j]-dx,xl[j],
                           rtol=1e-2)[1] / dx
    end

    # takes care of rounding error
    #can_opt.Pso[can_opt.Pso.>can_opt.Po]= minimum([can_opt.Po[can_opt.Pso.>can_opt.Po] can_opt.Ps[can_opt.Pso.>can_opt.Po]],dims=2)
    #can_opt.Pso[can_opt.Pso.>can_opt.Ps]= minimum([can_opt.Po[can_opt.Pso.>can_opt.Ps] can_opt.Ps[can_opt.Pso.>can_opt.Ps]],dims=2)

    return nothing
end
