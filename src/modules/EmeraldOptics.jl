module EmeraldOptics

using ..EmeraldConstants: AVOGADRO, H_PLANCK, LIGHT_SPEED


include("../../packages/EmeraldOptics.jl/src/photon.jl"       )
include("../../packages/EmeraldOptics.jl/src/transmittance.jl")

end
