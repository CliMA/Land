using Revise
using Land.Photosynthesis

const FT = Float32

met = meteo{Float32}()
l = Photosynthesis.leaf_params{Float32}()

mods = Photosynthesis.PhotoMods(
    fluorescence    = FlexasTolBerryFluorescence{FT}(),
    photosynthesis  = C3FvCBPhotoGs(),
    respiration     = RespirationCLM{FT}(),
    stomatal        = BallBerryStomata{FT}(),
    Jmax            = JmaxCLM{FT}(),
    Vmax            = VcmaxCLM{FT}(),
    MichaelisMenten = MM_CLM{FT}(),
    BoundaryLayer   = GentineLeafBoundary())

l.Kn = 2.44; l.α=0.2; l.ε=0.98; l.LMA=100e-3; l.RWC=80/100;l.psi_l=-1e6;l.psi_l50 = -1e6;l.ck=3;met.zscreen = 2.0;
l.height   = 1.0; met.zscreen  = 2.0;
l.dynamic_state = true

met.stab_type_stable = 2;
met.e_air = 1500;
met.T_air = 280;
Photosynthesis.LeafPhotosynthesis!(mods, l,met)