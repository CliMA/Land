# Photosynthesis.jl

Photosynthesis models for C3 and C4 photosynthesis.

## Usage
```julia
using Photosynthesis

mod_3 = C3CLM(Float32);
leaf  = Leaf{Float32}();
envir = PM.AirLayer{Float32}();
leaf_photo_from_pi(mod_3, leaf, envir);
```
