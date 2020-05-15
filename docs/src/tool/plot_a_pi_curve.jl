using LaTeXStrings
using PyPlot
using Land.Plant

# this function plot An-Pi curve
function plot_a_pi_curve(;
                         v25::FT = FT(80.0),
                         j25::FT = FT(135.0),
                      Γ_star::FT = FT(2.5),
                         tem::FT = FT(298.15),
                         par::FT = FT(1200.0),
                        p_O₂::FT = FT(21278.25),
                         r25::FT = FT(Inf)) where FT
    list_p,list_an,list_ag = Plant.get_a_pi_curve(
                                                  v25 = v25,
                                                  j25 = j25,
                                               Γ_star = Γ_star,
                                                  tem = tem,
                                                  par = par,
                                                 p_O₂ = p_O₂,
                                                  r25 = r25)

    # plot the data
    clf()
    plot(list_p, list_ag, color=:red  , label=L"$A_\mathrm{gross}$")
    plot(list_p, list_an, color=:green, label=L"$A_\mathrm{net}$"  )
    xlabel(L"$P_\mathrm{i}$ (Pa)"                                   , fontsize=16)
    ylabel(L"$A_\mathrm{net}$ ($\mathrm{\mu}$mol m$^{-2}$ s$^{-1}$)", fontsize=16)
    legend(loc="lower right")
end
