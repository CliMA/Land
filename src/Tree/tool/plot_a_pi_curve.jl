using LaTeXStrings
using PyPlot

# this function plot An-Pi curve
function plot_a_pi_curve(v25=80.0, j25=135.0, gamma=2.5, tem=298.15, par=1200.0, p_o2=21278.25, r25=false, unit="K")
    list_p,list_an,list_ag = get_a_pi_curve(v25, j25, gamma, tem, par, p_o2, r25, unit)

    # plot the data
    clf()
    plot(list_p, list_ag, color=:red  , label=L"$A_\mathrm{gross}$")
    plot(list_p, list_an, color=:green, label=L"$A_\mathrm{net}$"  )
    xlabel(L"$P_\mathrm{i}$ (Pa)"                                   , fontsize=16)
    ylabel(L"$A_\mathrm{net}$ ($\mathrm{\mu}$mol m$^{-2}$ s$^{-1}$)", fontsize=16)
    legend(loc="lower right")
end
