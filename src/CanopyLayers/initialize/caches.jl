###############################################################################
#
# Initialize CFCache
#
###############################################################################
"""
    create_cf_cache(FT, rt_dim::RTDimensions)

Create a [`CFCache`](@ref) type struct, given
- `FT` Floating number type
- `rt_dim` [`RTDimensions`](@ref) type struct
"""
function create_cf_cache(FT, rt_dim::RTDimensions)
    @unpack nAzi, nLayer, nPAR, nWL = rt_dim;

    abs_wave    = zeros(FT, nWL   );
    absfs_lidf  = zeros(FT, nAzi  );
    E_all       = zeros(FT, nWL   );
    E_iPAR      = zeros(FT, nPAR  );
    lPs         = zeros(FT, nLayer);
    kChlrel     = zeros(FT, nPAR  );
    PAR_diff    = zeros(FT, nPAR  );
    PAR_diffCab = zeros(FT, nPAR  );
    PAR_dir     = zeros(FT, nPAR  );
    PAR_dirCab  = zeros(FT, nPAR  );

    return CFCache{FT}(abs_wave    = abs_wave   ,
                       absfs_lidf  = absfs_lidf ,
                       E_all       = E_all      ,
                       E_iPAR      = E_iPAR     ,
                       lPs         = lPs        ,
                       kChlrel     = kChlrel    ,
                       PAR_diff    = PAR_diff   ,
                       PAR_diffCab = PAR_diffCab,
                       PAR_dir     = PAR_dir    ,
                       PAR_dirCab  = PAR_dirCab )
end








###############################################################################
#
# Initialize CGCache
#
###############################################################################
"""
    create_cg_cache(FT, rt_dim::RTDimensions)

Create a [`CGCache`](@ref) type struct, given
- `FT` Floating number type
- `rt_dim` [`RTDimensions`](@ref) type struct
"""
function create_cg_cache(FT, rt_dim::RTDimensions)
    @unpack nAzi, nIncl = rt_dim;

    _Co  = zeros(FT, nIncl);
    _Cs  = zeros(FT, nIncl);
    _So  = zeros(FT, nIncl);
    _Ss  = zeros(FT, nIncl);

    _1s  =  ones(FT, (1    , nAzi));
    _2d  = zeros(FT, (nIncl, nAzi));
    _cdo = zeros(FT, (nIncl, nAzi));
    _cds = zeros(FT, (nIncl, nAzi));

    return CGCache{FT}(_Co  = _Co ,
                       _Cs  = _Cs ,
                       _So  = _So ,
                       _Ss  = _Ss ,
                       _1s  = _1s ,
                       _2d  = _2d ,
                       _cdo = _cdo,
                       _cds = _cds)
end








###############################################################################
#
# Initialize SWCache
#
###############################################################################
"""
    create_sf_cache(FT, rt_dim::RTDimensions)

Create a [`SFCache`](@ref) type struct, given
- `FT` Floating number type
- `rt_dim` [`RTDimensions`](@ref) type struct
"""
function create_sf_cache(FT, rt_dim::RTDimensions)
    @unpack nAzi, nIncl, nLayer, nLevel, nWLE, nWLF = rt_dim;

    τ_dd            = zeros(FT, (nWLF,nLayer));
    Xdd             = zeros(FT, (nWLF,nLayer));
    Rdd             = zeros(FT, (nWLF,nLevel));
    Y               = zeros(FT, (nWLF,nLayer));
    U               = zeros(FT, (nWLF,nLevel));
    dnorm           = zeros(FT, nWLF);
    ρ_dd            = zeros(FT, (nWLF,nLayer));
    S⁻              = zeros(FT, (nWLF,nLayer));
    S⁺              = zeros(FT, (nWLF,nLayer));
    piLs            = zeros(FT, (nWLF,nLayer));
    piLd            = zeros(FT, (nWLF,nLayer));
    Fsmin           = zeros(FT, (nWLF,nLayer));
    Fsplu           = zeros(FT, (nWLF,nLayer));
    Fdmin           = zeros(FT, (nWLF,nLayer));
    Fdplu           = zeros(FT, (nWLF,nLayer));
    Femo            = zeros(FT, (nWLF,nLayer));
    M⁺              = zeros(FT, (nWLF,nWLE));
    M⁻              = zeros(FT, (nWLF,nWLE));
    M⁻_sun          = zeros(FT, nWLF);
    M⁺_sun          = zeros(FT, nWLF);
    wfEs            = zeros(FT, nWLF);
    sfEs            = zeros(FT, nWLF);
    sbEs            = zeros(FT, nWLF);
    M⁺⁻             = zeros(FT, nWLF);
    M⁺⁺             = zeros(FT, nWLF);
    M⁻⁺             = zeros(FT, nWLF);
    M⁻⁻             = zeros(FT, nWLF);
    sun_dwl_iWlE    = zeros(FT, nWLE);
    tmp_dwl_iWlE    = zeros(FT, nWLE);
    ϕ_cosΘ          = zeros(FT, (nIncl,nAzi));
    ϕ_cosΘ_lidf     = zeros(FT, nAzi);
    vfEplu_shade    = zeros(FT, nWLF);
    vbEmin_shade    = zeros(FT, nWLF);
    vfEplu_sun      = zeros(FT, nWLF);
    vbEmin_sun      = zeros(FT, nWLF);
    sigfEmin_shade  = zeros(FT, nWLF);
    sigbEmin_shade  = zeros(FT, nWLF);
    sigfEmin_sun    = zeros(FT, nWLF);
    sigbEmin_sun    = zeros(FT, nWLF);
    sigfEplu_shade  = zeros(FT, nWLF);
    sigbEplu_shade  = zeros(FT, nWLF);
    sigfEplu_sun    = zeros(FT, nWLF);
    sigbEplu_sun    = zeros(FT, nWLF);
    zeroB           = zeros(FT, nWLF);
    F⁻              = zeros(FT, (nWLF,nLevel));
    F⁺              = zeros(FT, (nWLF,nLevel));
    net_diffuse     = zeros(FT, (nWLF,nLayer));
    tmp_1d_nWlF     = zeros(FT, nWLF);
    tmp_1d_nLayer   = zeros(FT, nLayer);
    tmp_2d_nWlF_nLayer   = zeros(FT, (nWLF,nLayer));
    tmp_2d_nWlF_nLayer_2 = zeros(FT, (nWLF,nLayer));

    return SFCache{FT}(τ_dd                 = τ_dd                ,
                       Xdd                  = Xdd                 ,
                       Rdd                  = Rdd                 ,
                       Y                    = Y                   ,
                       U                    = U                   ,
                       dnorm                = dnorm               ,
                       ρ_dd                 = ρ_dd                ,
                       S⁻                   = S⁻                  ,
                       S⁺                   = S⁺                  ,
                       piLs                 = piLs                ,
                       piLd                 = piLd                ,
                       Fsmin                = Fsmin               ,
                       Fsplu                = Fsplu               ,
                       Fdmin                = Fdmin               ,
                       Fdplu                = Fdplu               ,
                       Femo                 = Femo                ,
                       M⁺                   = M⁺                  ,
                       M⁻                   = M⁻                  ,
                       M⁻_sun               = M⁻_sun              ,
                       M⁺_sun               = M⁺_sun              ,
                       wfEs                 = wfEs                ,
                       sfEs                 = sfEs                ,
                       sbEs                 = sbEs                ,
                       M⁺⁻                  = M⁺⁻                 ,
                       M⁺⁺                  = M⁺⁺                 ,
                       M⁻⁺                  = M⁻⁺                 ,
                       M⁻⁻                  = M⁻⁻                 ,
                       sun_dwl_iWlE         = sun_dwl_iWlE        ,
                       tmp_dwl_iWlE         = tmp_dwl_iWlE        ,
                       ϕ_cosΘ               = ϕ_cosΘ              ,
                       ϕ_cosΘ_lidf          = ϕ_cosΘ_lidf         ,
                       vfEplu_shade         = vfEplu_shade        ,
                       vbEmin_shade         = vbEmin_shade        ,
                       vfEplu_sun           = vfEplu_sun          ,
                       vbEmin_sun           = vbEmin_sun          ,
                       sigfEmin_shade       = sigfEmin_shade      ,
                       sigbEmin_shade       = sigbEmin_shade      ,
                       sigfEmin_sun         = sigfEmin_sun        ,
                       sigbEmin_sun         = sigbEmin_sun        ,
                       sigfEplu_shade       = sigfEplu_shade      ,
                       sigbEplu_shade       = sigbEplu_shade      ,
                       sigfEplu_sun         = sigfEplu_sun        ,
                       sigbEplu_sun         = sigbEplu_sun        ,
                       zeroB                = zeroB               ,
                       F⁻                   = F⁻                  ,
                       F⁺                   = F⁺                  ,
                       net_diffuse          = net_diffuse         ,
                       tmp_1d_nWlF          = tmp_1d_nWlF         ,
                       tmp_1d_nLayer        = tmp_1d_nLayer       ,
                       tmp_2d_nWlF_nLayer   = tmp_2d_nWlF_nLayer  ,
                       tmp_2d_nWlF_nLayer_2 = tmp_2d_nWlF_nLayer_2)
end








###############################################################################
#
# Initialize SWCache
#
###############################################################################
"""
    create_sw_cache(FT, rt_dim::RTDimensions)

Create a [`CGCache`](@ref) type struct, given
- `FT` Floating number type
- `rt_dim` [`RTDimensions`](@ref) type struct
"""
function create_sw_cache(FT, rt_dim::RTDimensions)
    @unpack nLayer, nWL = rt_dim;

    dnorm  = zeros(FT, nWL);
    piLo   = zeros(FT, nWL);
    piLoc  = zeros(FT, nWL);
    piLos  = zeros(FT, nWL);

    piLoc2 = zeros(FT, (nWL,nLayer));
    ρ_dd   = zeros(FT, (nWL,nLayer));
    ρ_sd   = zeros(FT, (nWL,nLayer));
    τ_dd   = zeros(FT, (nWL,nLayer));
    τ_sd   = zeros(FT, (nWL,nLayer));

    return SWCache{FT}(dnorm  = dnorm ,
                       piLo   = piLo  ,
                       piLoc  = piLoc ,
                       piLos  = piLos ,
                       piLoc2 = piLoc2,
                       ρ_dd   = ρ_dd  ,
                       ρ_sd   = ρ_sd  ,
                       τ_dd   = τ_dd  ,
                       τ_sd   = τ_sd  )
end








###############################################################################
#
# Initialize RTCache
#
###############################################################################
"""
    create_rt_cache(FT, rt_dim::RTDimensions)

Create an [`RTCache`](@ref), given
- `FT` Floating number type
- `rt_dim` [`RTDimensions`](@ref) type struct
"""
function create_rt_cache(FT, rt_dim::RTDimensions)
    cf_con = create_cf_cache(FT, rt_dim);
    cg_con = create_cg_cache(FT, rt_dim);
    sf_con = create_sf_cache(FT, rt_dim);
    sw_con = create_sw_cache(FT, rt_dim);

    return RTCache{FT}(cf_con = cf_con,
                       cg_con = cg_con,
                       sf_con = sf_con,
                       sw_con = sw_con)
end
