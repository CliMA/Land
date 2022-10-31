# StomtaModels.jl

## Use StomataModels
```
julia> using Photosynthesis
julia> using StomataModels
julia> envir  = AirLayer{FT}();
julia> ps_3   = C3CLM(FT);
julia> leaves = Leaves{FT}(n_leaf=2);
julia> sm     = OSMWang{FT}();
julia>
julia> gas_exchange!(ps_3, leaves, envir, sm);
```
