# Photosynthesis.jl

Photosynthesis models for C3 and C4 photosynthesis.

## Usage
```julia
using Photosynthesis

mod_3 = C3CLM(Float32);
leaf  = Leaf{Float32}();
envir = AirLayer{Float32}();
leaf_photosynthesis!(mod_3, leaf, envir, GCO₂Mode());
```
