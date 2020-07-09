function Yujie111UpdateLegacy(node::Yujie111{FT}, f_sl, f_sh, r_sl) where {FT<:AbstractFloat}
    # 0.0 sum the flow for rhiz. root, and stem
    flow = f_sl + f_sh

    # 1.1 assign the soil parameters
    p_ssat  = node.p_ssat
    c_ssat  = node.c_ssat
    b_ssat  = node.b_ssat
    k_rhiz  = node.k_rhiz
    p_soil  = node.p_soil
    tension = node.p_soil
    dp      = 0.0
    # 1.2 iterate through each shell to get dp
    for i = 0:9
        if tension>p_ssat
            shell_t = c_ssat * (tension/p_ssat)^(-1.0/b_ssat)
        else
            shell_t = c_ssat
        end
        shell_f  = (shell_t/c_ssat) ^ (2.0*b_ssat+3)
        shell_k  = k_rhiz * shell_f *
                   log(10.0) / log((10.0-0.9*i)/(10.0-0.9*(i+1)))
        dp      += flow / shell_k
        tension  = p_soil + dp
    end
    p_rhiz  = tension
    tension = p_rhiz


    # 2.1 assign the root parameters
    vc      = node.vc_root
    k_root  = node.k_root
    h       = node.h_root
    legacy  = node.l_root
    dp      = 0.0
    # 2.2 iterate through each layer to get root dp including gravity
    for i in 1:20
        if tension>legacy[i,1]
            layer_f = xylem_k_ratio(vc,-tension)
            node.l_root[i,1] = tension
            node.l_root[i,2] = layer_f
        else
            layer_f = legacy[i,2]
        end
        layer_k  = k_root * layer_f * 20.0
        dp      += flow/layer_k + 998.0*9.8*h*0.05*1E-6
        tension  = p_rhiz + dp
    end
    p_root  = tension
    tension = p_root


    # 3.1 assign the stem parameters
    vc      = node.vc_stem
    k_stem  = node.k_stem
    h       = node.h_stem
    legacy  = node.l_stem
    dp      = 0.0
    # 3.2 iterate through each layer to get stem dp including gravity
    for i in 1:20
        if tension>legacy[i,1]
            layer_f = xylem_k_ratio(vc,-tension)
            node.l_stem[i,1] = tension
            node.l_stem[i,2] = layer_f
        else
            layer_f = legacy[i,2]
        end
        layer_k  = k_stem * layer_f * 20.0
        dp      += flow/layer_k + 998.0*9.8*h*0.05*1E-6
        tension  = p_root + dp
    end
    p_stem  = tension
    tension = p_stem


    # 4.1 assign the sunlit leaf parameters
    vc      = node.vc_leaf
    k_leaf  = node.k_leaf * r_sl
    h       = node.h_leaf
    legacy  = node.l_leaf
    dp      = 0.0
    # 4.2 iterate through each sunlit layer to get leaf dp including gravity
    for i in 1:20
        if tension>legacy[i,1]
            layer_f = xylem_k_ratio(vc,-tension)
            node.l_leaf[i,1] = tension
            node.l_leaf[i,2] = layer_f
        else
            layer_f = legacy[i,2]
        end
        layer_k  = k_leaf * layer_f * 20.0
        dp      += f_sl/layer_k + 998.0*9.8*h*0.05*1E-6
        tension  = p_stem + dp
    end


    # 5.1 assign the shade leaf parameters
    tension = p_stem
    k_leaf  = node.k_leaf * (1.0-r_sl)
    dp      = 0.0
    # 5.2 iterate through each shade layer to get leaf dp including gravity
    for i in 1:20
        if tension>legacy[i,1]
            layer_f = xylem_k_ratio(vc,-tension)
            node.l_leaf[i,1] = tension
            node.l_leaf[i,2] = layer_f
        else
            layer_f = legacy[i,2]
        end
        layer_k  = k_leaf * layer_f * 20.0
        dp      += f_sh/layer_k + 998.0*9.8*h*0.05*1E-6
        tension  = p_stem + dp
    end
end
