"""
    sif_fluxes!(leaf_array::Array{LeafBios{FT},1}, can_opt::CanopyOpticals{FT}, can_rad::CanopyRads{FT}, can::Canopy4RT{FT}, soil_opt::SoilOpticals{FT}, wl_set::WaveLengths{FT}) where {FT<:AbstractFloat}

Computes 2-stream diffusive radiation transport for SIF radiation (calls `compute_diffusive_S` internally).
Layer reflectance and transmission is computed from SW optical properties, layer sources from absorbed light and SIF efficiencies. Boundary conditions are zero SIF incoming from atmosphere or soil.
- `leaf_array` An array of [`LeafBioArray`](@ref) type struct (i.e. leaf optical properties can change with canopy height)
- `can_opt` A [`CanopyOptiArray`](@ref) struct for providing optical layer properties
- `can_rad` A [`CanopyRadiation`](@ref) struct
- `can` A [`Canopy4RT`](@ref) type struct for providing LAI and nlayers and clumping
- `soil_opt` A [`SoilOpti`](@ref) type struct for soil optical properties
- `wl_set` An [`WLParaSetArray`](@ref) type struct
"""
function sif_fluxes!(
            leaf_array::Array{LeafBios{FT},1},
            can_opt::CanopyOpticals{FT},
            can_rad::CanopyRads{FT},
            can::Canopy4RT{FT},
            soil_opt::SoilOpticals{FT},
            wl_set::WaveLengths{FT}
) where {FT<:AbstractFloat}
    @unpack fs,fo, cosΘ_l, absfs, fsfo, absfsfo, cos2Θ_l, Ps, Po, Pso, a,  sigb,vb,vf = can_opt
    @unpack E_down, E_up,ϕ_shade,ϕ_sun = can_rad
    @unpack Ω,nlayers, LAI = can
    @unpack dwl, Iwle, Iwlf = wl_set
    rsoil = soil_opt.albedo_SW[Iwlf]

    iLAI  = Ω*LAI/nlayers;

    τ_dd  = 1 .-a[Iwlf,:]*iLAI;
    ρ_dd  = sigb[Iwlf,:]*iLAI;

    dim   = size(leaf_array[1].Mb)
    # Compute mid layer Ps,Po,Pso
    #Qso   = (Pso[1:nlayers] + Pso[2:nlayers+1])/2;
    #Qs      = (Ps[1:nlayers] + Ps[2:nlayers+1])/2;
    #Qo      =(Po[1:nlayers] + Po[2:nlayers+1])/2;
    Qso   = Pso[1:nlayers]
    Qs    = Ps[1:nlayers]
    Qo    = Po[1:nlayers]
    # Allocate arrays for output:
    di    = size(leaf_array[1].Mb)
    S⁻    = zeros(FT,dim[1],nlayers )
    S⁺    = similar(S⁻)
    piLs  = similar(S⁻)
    piLd  = similar(S⁻)
    Fsmin = similar(S⁻)
    Fsplu = similar(S⁻)
    Fdmin = similar(S⁻)
    Fdplu = similar(S⁻)
    Femo  = similar(S⁻);

    # fs is different here than in RTMo (in RTMo, it is abs(fs))!

    # 2. Calculation of reflectance
    # 2.1  reflectance, transmittance factors in a thin layer the following are vectors with length [nl,nwl]
    sun = can_opt.Es_[Iwle,1];
    @inbounds for i=1:nlayers
        if length(leaf_array)>1
            Mb = leaf_array[i].Mb
            Mf = leaf_array[i].Mf
        else
            Mb = leaf_array[1].Mb
            Mf = leaf_array[1].Mf
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
        sunCos    = mean((ϕ_sun[:,:,i].*cosΘ_l)'*can.lidf)
        shadeCos  = mean((ϕ_shade[i]*cosΘ_l)'*can.lidf)
        sunCos2   = mean((ϕ_sun[:,:,i].*cos2Θ_l)'*can.lidf)
        shadeCos2 = mean((ϕ_shade[i]*cos2Θ_l)'*can.lidf)
        sunLidf   = mean(ϕ_sun[:,:,i]'*can.lidf)
        shadeLidf = mean(ϕ_shade[i]'*can.lidf)

        wfEs = mean((ϕ_sun[:,:,i].*absfsfo)'*can.lidf)  * M⁺_sun + mean((ϕ_sun[:,:,i].*fsfo)'*can.lidf)  *M⁻_sun

        a1   = mean((ϕ_sun[:,:,i].*absfs)'*can.lidf)* M⁺_sun;
        a2   = mean((ϕ_sun[:,:,i].*fs.*cosΘ_l)'*can.lidf)*M⁻_sun
        sfEs = a1 - a2
        sbEs = a1 + a2

        a1 = mean((ϕ_shade[i]*abs.(fo))'*can.lidf) ;
        a2 = mean((ϕ_shade[i]*fo.*cosΘ_l)'*can.lidf)
        vfEplu_shade =  a1  * M⁺⁺ + a2 *M⁻⁺
        vbEmin_shade =  a1 * M⁺⁻ + a2 *M⁻⁻

        a1 = mean((ϕ_sun[:,:,i].*abs.(fo))'*can.lidf);
        a2 = mean((ϕ_sun[:,:,i].*fo.*cosΘ_l)'*can.lidf)
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
    F⁻,F⁺,net_diffuse = diffusive_S(τ_dd, ρ_dd,S⁻, S⁺,zeroB, zeroB, rsoil)
    # Save in output structures!
    can_rad.SIF_obs_sunlit[:] = iLAI/FT(pi)*Qso'*piLs';                                               # direct Sunlit leaves
    can_rad.SIF_obs_shaded[:] = iLAI/FT(pi)*(Qo[1:nlayers]-Qso[1:nlayers])'*piLd';                    # direct shaded leaves

    # SIF scattered internally
    #can_rad.SIF_obs_scattered[:]     = iLAI/FT(pi)*(Qo[1:nlayers]'*(vb[Iwlf,:].*F⁻[:,1:nlayers] + vf[Iwlf,:].*F⁺[:,1:nlayers])');
    can_rad.SIF_obs_scattered[:] = iLAI/FT(pi)*(Qo[1:nlayers]'*(vb[Iwlf,:].*F⁻[:,1:nlayers] + vf[Iwlf,:].*F⁺[:,1:nlayers])');
    can_rad.SIF_obs_soil[:]      = (rsoil .* F⁻[:,end] * Po[end])/FT(pi); #Soil contribution

    can_rad.SIF_hemi[:] = F⁺[:,1];
    can_rad.SIF_obs[:]  = can_rad.SIF_obs_sunlit[:]+can_rad.SIF_obs_shaded[:]+can_rad.SIF_obs_scattered[:] +can_rad.SIF_obs_soil[:];
    can_rad.SIF_sum[:]  = sum(S⁻+S⁺, dims=2)
    #return  F⁻,F⁺,S⁻,S⁺, piLs, piLd

    return nothing
end
