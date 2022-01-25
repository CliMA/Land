# Photosynthesis.jl

Photosynthesis models for C3 and C4 photosynthesis.

## Install
```julia
using Pkg;
Pkg.add("Photosynthesis");
```


## Examples
```julia
leaf   = Leaf{Float32}("C3");
air    = AirLayer{Float32}();
p_mode = PCO₂Mode();
g_mode = PCO₂Mode();
leaf_photosynthesis!(leaf, air, p_mode);
leaf_photosynthesis!(leaf, air, g_mode);
```
