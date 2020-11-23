# PlantHydraulics.jl
Plant hydraulics package using numerical methods.

## Usage
```julia
using PlantHydraulics

leaf   = LeafHydraulics{Float32}();
critical_flow = critical_flow(leaf);
```
