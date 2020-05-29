# this function plot An-Pi curve
function plot_a_par_curve(;
                          v25::FT = FT(80.0),
                          j25::FT = FT(135.0),
                       Γ_star::FT = FT(2.5),
                          gsc::FT = FT(0.1),
                          p_a::FT = FT(40.0),
                          tem::FT = FT(298.15),
                        p_atm::FT = FT(101325.0),
                         p_O₂::FT = FT(21278.25),
                          r25::FT = FT(Inf)) where {FT}
    list_par,list_an,list_ag,list_pi = get_a_par_curve(
                                                       v25 = v25,
                                                       j25 = j25,
                                                    Γ_star = Γ_star,
                                                       gsc = gsc,
                                                       p_a = p_a,
                                                       tem = tem,
                                                     p_atm = p_atm,
                                                      p_O₂ = p_O₂,
                                                       r25 = r25)

    # plot the data
    clf()
    plot(list_par, list_ag, color=:red  , label=L"$A_\mathrm{gross}$")
    plot(list_par, list_an, color=:green, label=L"$A_\mathrm{net}$"  )
    xlabel(L"PAR ($\mathrm{\mu}$mol m$^{-2}$ s$^{-1}$)"             , fontsize=16)
    ylabel(L"$A_\mathrm{net}$ ($\mathrm{\mu}$mol m$^{-2}$ s$^{-1}$)", fontsize=16)
    legend(loc="center right")

    twinx()
    plot(list_par, list_pi, color=:blue , label=L"$P_\mathrm{i}$"    )
    ylabel(L"$P_\mathrm{i}$ (Pa)", fontsize=16, color=:blue)
end
