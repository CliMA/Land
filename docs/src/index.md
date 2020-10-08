# WaterPhysics.jl

Temperature dependencies of water and other trace molecule physical properties. See [`WaterPhysics.jl`](https://github.com/Yujie-W/WaterPhysics.jl) for the source code.





## Installation
```julia
julia> using Pkg;
julia> Pkg.add("WaterPhysics");
```




## Example
```julia
using WaterPhysics

# calculate the relative viscosity of liquid water at 300 K
f_vis = relative_viscosity(300.0);
```
