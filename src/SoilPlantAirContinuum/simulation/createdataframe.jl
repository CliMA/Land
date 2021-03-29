###############################################################################
#
# Create DataFrame to save simulation output
#
###############################################################################
"""
    function create_dataframe(FT, weather::DataFrame)

Create a data frame to store simulation output, given
- `FT` Floating number type
- `weather` Weather profile in a growing season
"""
function create_dataframe(
            FT,
            weather::DataFrame
)
    df    = DataFrame();

    # climatic info
    df[!, "Time"  ]  = (weather).Day + (weather).Hour / FT(24);
    df[!, "T_air" ] .= FT(0);
    df[!, "D_air" ] .= FT(0);
    df[!, "Wind"  ] .= FT(0);
    df[!, "Rain"  ] .= FT(0);
    df[!, "Ca"    ] .= FT(0);

    # whole tree level
    df[!, "SWC"   ] .= FT(0);
    df[!, "P_soil"] .= FT(0);
    df[!, "H_sun" ] .= FT(0);
    df[!, "A_net" ] .= FT(0);
    df[!, "E_crit"] .= FT(0);

    # sunlit layers
    df[!, "LAI_sl"] .= FT(0);
    df[!, "PAR_sl"] .= FT(0);
    df[!, "RAD_sl"] .= FT(0);
    df[!, "E_sl"  ] .= FT(0);
    df[!, "P_sl"  ] .= FT(0);
    df[!, "An_sl" ] .= FT(0);
    df[!, "Ag_sl" ] .= FT(0);
    df[!, "C_sl"  ] .= FT(0);
    df[!, "G_sl"  ] .= FT(0);
    df[!, "T_sl"  ] .= FT(0);

    # shaded layers
    df[!, "LAI_sh"] .= FT(0);
    df[!, "PAR_sh"] .= FT(0);
    df[!, "RAD_sh"] .= FT(0);
    df[!, "E_sh"  ] .= FT(0);
    df[!, "P_sh"  ] .= FT(0);
    df[!, "An_sh" ] .= FT(0);
    df[!, "Ag_sh" ] .= FT(0);
    df[!, "C_sh"  ] .= FT(0);
    df[!, "G_sh"  ] .= FT(0);
    df[!, "T_sh"  ] .= FT(0);

    return df
end
