function Yujie111GetLeafPartition(node::Yujie111{FT}, zenith, r_all) where {FT<:AbstractFloat}
    # 1. assume 50%-50% PAR and NIR, 80%-20% beam-diffuse
    par_tot = r_all * 0.5 / 0.235
    q_ob    = par_tot * 0.8
    q_od    = par_tot * 0.2
    k_d     = 0.7
    alpar   = 0.8
    alnir   = 0.2


    # 2. calculate the LAI
    shape   = 1.0
    shapa   = 2.0 # hemi-sphere
    lai_t   = node.laba / node.gaba


    # 3. calculate the mean sunlit layer PAR
    # For vertical leaves k_be = 2.0 * tand(zenith) / pi
    k_be    = sqrt( shape^2+tand(zenith)^2 ) /
              ( shape+1.774*(shape+1.182)^(-0.733) )
    q_dp    = q_od * exp( -1.0*sqrt(alpar)*k_d *lai_t )
    q_dn    = q_od * exp( -1.0*sqrt(alnir)*k_d *lai_t )
    q_btp   = q_ob * exp( -1.0*sqrt(alpar)*k_be*lai_t )
    q_btn   = q_ob * exp( -1.0*sqrt(alnir)*k_be*lai_t )
    q_bp    = q_ob * exp( -1.0*k_be*lai_t             )
    q_bn    = q_ob * exp( -1.0*k_be*lai_t             )
    q_scp   = q_btp - q_bp
    q_scn   = q_btn - q_bn


    # 4. mean par in sunlit and shade layers
    q_dm    = q_od * ( 1.0-exp(-1.0*sqrt(alpar)*k_d) ) / ( sqrt(alpar)*k_d*lai_t )
    r_dm    = q_od * ( 1.0-exp(-1.0*sqrt(alnir)*k_d) ) / ( sqrt(alnir)*k_d*lai_t )
    q_shm   = q_dm + q_scp/shapa
    r_shm   = alnir * ( r_dm + q_scn/shapa )
    q_slm   = k_be*q_ob + q_dm + q_scp/shapa
    r_slm   = alnir * ( k_be*q_ob + r_dm + q_scn/shapa )


    # 5. calculate lasi_sl and lai_sh
    lai_sl = ( 1.0-exp(-1.0*k_be*lai_t) ) / k_be
    lai_sh = lai_t - lai_sl
    ratio = lai_sl / lai_t


    # 6. calculate the radiation based on PAR and NIR
    ep_sl  = q_slm * alpar * 0.235
    ep_sh  = q_shm * alpar * 0.235
    en_sl  = r_all * 0.5 * r_slm / (r_slm+r_shm)
    en_sh  = r_all * 0.5 * r_shm / (r_slm+r_shm)
    # e_sl and e_sh are layer-based radiation
    e_sl   = ep_sl*lai_sl + en_sl
    e_sh   = ep_sh*lai_sh + en_sh

    # 7. return the information
    return [ratio, q_slm, e_sl, 1.0-ratio, q_shm, e_sh]
end
