# # Predefined parameter sets

## load packages
using Land.Photosynthesis
using PlotPlants
FT = Float32;
#------------------------------------------------------------------------------




# ## Jmax
_td_1 = JmaxTDBernacchi(FT);
_td_2 = JmaxTDCLM(FT);
_td_3 = JmaxTDLeuning(FT);
_ts   = collect(FT, 273:1:323);
_jm_1 = photo_TD_from_val.([_td_1], FT(100), _ts);
_jm_2 = photo_TD_from_val.([_td_2], FT(100), _ts);
_jm_3 = photo_TD_from_val.([_td_3], FT(100), _ts);

_fig,_axes = create_canvas("Jmax");
_ax1 = _axes[1];
_ax1.plot(_ts .- 273.15, _jm_1, "k-", label="Bernacchi");
_ax1.plot(_ts .- 273.15, _jm_2, "k:", label="CLM");
_ax1.plot(_ts .- 273.15, _jm_3, "k--", label="Leuning");
_ax1.set_xlabel("Leaf temperature (°C)");
_ax1.set_ylabel("Jcmax (μmol m⁻² s⁻¹)");
_ax1.legend(loc="upper left");
_fig
#------------------------------------------------------------------------------




# ## Kc
_td_1 = KcTDBernacchi(FT);
_td_2 = KcTDCLM(FT);
_ts   = collect(FT, 273:1:323);
_kc_1 = photo_TD_from_set.([_td_1], _ts);
_kc_2 = photo_TD_from_set.([_td_2], _ts);

_fig,_axes = create_canvas("Kc");
_ax1 = _axes[1];
_ax1.plot(_ts .- 273.15, _kc_1, "k-", label="Bernacchi");
_ax1.plot(_ts .- 273.15, _kc_2, "k:", label="CLM");
_ax1.set_xlabel("Leaf temperature (°C)");
_ax1.set_ylabel("Kc (Pa)");
_ax1.legend(loc="upper left");
_fig
#------------------------------------------------------------------------------




# ## Ko
_td_1 = KoTDBernacchi(FT);
_td_2 = KoTDCLM(FT);
_ts   = collect(FT, 273:1:323);
_ko_1 = photo_TD_from_set.([_td_1], _ts);
_ko_2 = photo_TD_from_set.([_td_2], _ts);

_fig,_axes = create_canvas("Ko");
_ax1 = _axes[1];
_ax1.plot(_ts .- 273.15, _ko_1, "k-", label="Bernacchi");
_ax1.plot(_ts .- 273.15, _ko_2, "k:", label="CLM");
_ax1.set_xlabel("Leaf temperature (°C)");
_ax1.set_ylabel("Ko (Pa)");
_ax1.legend(loc="upper left");
_fig
#------------------------------------------------------------------------------




# ## Kpep
_td_1 = KpepTDBoyd(FT);
_td_2 = KpepTDCLM(FT);
_ts   = collect(FT, 273:1:323);
_kp_1 = photo_TD_from_set.([_td_1], _ts);
_kp_2 = photo_TD_from_set.([_td_2], _ts);

_fig,_axes = create_canvas("Kpep");
_ax1 = _axes[1];
_ax1.plot(_ts .- 273.15, _kp_1, "k-", label="Boyd");
_ax1.plot(_ts .- 273.15, _kp_2, "k:", label="CLM");
_ax1.set_xlabel("Leaf temperature (°C)");
_ax1.set_ylabel("Kpep (Pa)");
_ax1.legend(loc="upper left");
_fig
#------------------------------------------------------------------------------




# ## Respiration
_td_1 = RespirationTDBernacchi(FT);
_td_2 = RespirationTDCLM(FT);
_td_3 = Q10TDAngiosperm(FT);
_td_4 = Q10TDGymnosperm(FT);
_ts   = collect(FT, 273:1:323);
_rd_1 = photo_TD_from_val.([_td_1], FT(1), _ts);
_rd_2 = photo_TD_from_val.([_td_2], FT(1), _ts);
_rd_3 = photo_TD_from_val.([_td_3], FT(1), _ts);
_rd_4 = photo_TD_from_val.([_td_4], FT(1), _ts);

_fig,_axes = create_canvas("Respiration");
_ax1 = _axes[1];
_ax1.plot(_ts .- 273.15, _rd_1, "k-", label="Bernacchi");
_ax1.plot(_ts .- 273.15, _rd_2, "k:", label="CLM");
_ax1.plot(_ts .- 273.15, _rd_3, "r-", label="Q10 Angiosperm");
_ax1.plot(_ts .- 273.15, _rd_4, "r:", label="Q10 Gymnosperm");
_ax1.set_xlabel("Leaf temperature (°C)");
_ax1.set_ylabel("Respiration (μmol m⁻² s⁻¹)");
_ax1.legend(loc="upper left");
_fig
#------------------------------------------------------------------------------




# ## Vcmax
_td_1 = VcmaxTDBernacchi(FT);
_td_2 = VcmaxTDCLM(FT);
_td_3 = VcmaxTDLeuning(FT);
_ts   = collect(FT, 273:1:323);
_vc_1 = photo_TD_from_val.([_td_1], FT(100), _ts);
_vc_2 = photo_TD_from_val.([_td_2], FT(100), _ts);
_vc_3 = photo_TD_from_val.([_td_3], FT(100), _ts);

_fig,_axes = create_canvas("Vcmax");
_ax1 = _axes[1];
_ax1.plot(_ts .- 273.15, _vc_1, "k-", label="Bernacchi");
_ax1.plot(_ts .- 273.15, _vc_2, "k:", label="CLM");
_ax1.plot(_ts .- 273.15, _vc_3, "k--", label="Leuning");
_ax1.set_xlabel("Leaf temperature (°C)");
_ax1.set_ylabel("Vcmax (μmol m⁻² s⁻¹)");
_ax1.legend(loc="upper left");
_fig
#------------------------------------------------------------------------------




# ## Vomax
_td_1 = VomaxTDBernacchi(FT);
_ts   = collect(FT, 273:1:323);
_vo_1 = photo_TD_from_val.([_td_1], FT(100), _ts);

_fig,_axes = create_canvas("Vomax");
_ax1 = _axes[1];
_ax1.plot(_ts .- 273.15, _vo_1, "k-", label="Bernacchi");
_ax1.set_xlabel("Leaf temperature (°C)");
_ax1.set_ylabel("Vomax (μmol m⁻² s⁻¹)");
_ax1.legend(loc="upper left");
_fig
#------------------------------------------------------------------------------




# ## Vpmax
_td_1 = VpmaxTDBoyd(FT);
_ts   = collect(FT, 273:1:323);
_vp_1 = photo_TD_from_val.([_td_1], FT(100), _ts);

_fig,_axes = create_canvas("Vpmax");
_ax1 = _axes[1];
_ax1.plot(_ts .- 273.15, _vp_1, "k-", label="Boyd");
_ax1.set_xlabel("Leaf temperature (°C)");
_ax1.set_ylabel("Vpmax (μmol m⁻² s⁻¹)");
_ax1.legend(loc="upper left");
_fig
#------------------------------------------------------------------------------




# ## Γ*
_td_1 = ΓStarTDBernacchi(FT);
_td_2 = ΓStarTDCLM(FT);
_ts   = collect(FT, 273:1:323);
_Γs_1 = photo_TD_from_set.([_td_1], _ts);
_Γs_2 = photo_TD_from_set.([_td_2], _ts);

_fig,_axes = create_canvas("Γ*");
_ax1 = _axes[1];
_ax1.plot(_ts .- 273.15, _Γs_1, "k-", label="Bernacchi");
_ax1.plot(_ts .- 273.15, _Γs_2, "k:", label="CLM");
_ax1.set_xlabel("Leaf temperature (°C)");
_ax1.set_ylabel("Γ* (Pa)");
_ax1.legend(loc="upper left");
_fig
#------------------------------------------------------------------------------
