###############################################################################
#
# Big leaf canopy model
#
###############################################################################
const _K_D         = 0.7
const _ALPAR       = 0.8
const _ALNIR       = 0.2
const _ALPAR_SQ    = sqrt(_ALPAR)
const _ALNIR_SQ    = sqrt(_ALNIR)
const _ALPAR_SQ_KD = _ALPAR_SQ * _K_D
const _ALNIR_SQ_KD = _ALNIR_SQ * _K_D




"""
    big_leaf_partition!(partition::Array{FT,1}, zenith::FT, lai::FT, r_all::FT)

Partition the big-leaf canopy into sunlit and shaded layers, given
- `partition` Container for partition
- `zenith` Zenith angle in degree
- `lai` Leaf area index
- `r_all` Total radiation in `[W m⁻²]`
"""
function big_leaf_partition!(
            partition::Array{FT,1},
            zenith::FT,
            lai::FT,
            r_all::FT
            ) where {FT <:AbstractFloat}
    # 1. assume 50%-50% PAR and NIR, 80%-20% beam-diffuse
    par_tot = r_all / 2 / FT(0.235);
    q_ob    = par_tot * FT(0.8);
    q_od    = par_tot * FT(0.2);

    # 2. calculate the LAI
    shape   = FT(1.0);
    shapa   = FT(2.0);

    # 3. calculate the mean sunlit layer PAR
    # For vertical leaves k_be = 2.0 * tand(zenith) / pi
    k_be    = sqrt( shape^2+tand(zenith)^2 ) /
              ( shape + FT(1.774) * (shape+FT(1.182))^FT(-0.733) );
    q_btp   = q_ob * exp( -FT(_ALPAR_SQ) * k_be * lai );
    q_btn   = q_ob * exp( -FT(_ALNIR_SQ) * k_be * lai );
    q_b     = q_ob * exp( -k_be*lai                   );
    q_scp   = q_btp - q_b;
    q_scn   = q_btn - q_b;

    # 4. mean par in sunlit and shade layers
    q_dm    = q_od * ( 1 - exp(-FT(_ALPAR_SQ_KD)) ) / ( FT(_ALPAR_SQ_KD) * lai );
    r_dm    = q_od * ( 1 - exp(-FT(_ALNIR_SQ_KD)) ) / ( FT(_ALNIR_SQ_KD) * lai );
    q_shm   = q_dm + q_scp/shapa;
    r_shm   = _ALNIR * ( r_dm + q_scn/shapa );
    q_slm   = k_be*q_ob + q_dm + q_scp/shapa;
    r_slm   = _ALNIR * ( k_be*q_ob + r_dm + q_scn/shapa );

    # 5. calculate lasi_sl and lai_sh
    lai_sl = ( 1 - exp(-k_be*lai) ) / k_be
    lai_sh = lai - lai_sl
    ratio = lai_sl / lai

    # 6. calculate the radiation based on PAR and NIR
    ep_sl  = q_slm * _ALPAR * FT(0.235);
    ep_sh  = q_shm * _ALPAR * FT(0.235);
    en_sl  = r_all / 2 * r_slm / (r_slm+r_shm);
    en_sh  = r_all / 2 * r_shm / (r_slm+r_shm);
    # e_sl and e_sh are layer-based radiation
    e_sl   = ep_sl*lai_sl + en_sl;
    e_sh   = ep_sh*lai_sh + en_sh;

    # 7. update the information
    partition[1] = ratio  ;
    partition[2] = q_slm  ;
    partition[3] = e_sl   ;
    partition[4] = 1-ratio;
    partition[5] = q_shm  ;
    partition[6] = e_sh   ;

    return nothing
end
