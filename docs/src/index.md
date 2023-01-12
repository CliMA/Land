# CliMA Land

## Install CliMA Land dev (v0.2)

```julia
julia> using Pkg
julia> Pkg.add(PackageSpec(url="https://github.com/CliMA/Land.git"))
```

## Install CliMA Land v0.1

```julia
julia> using Pkg
julia> Pkg.add(PackageSpec(url="https://github.com/CliMA/Land.git", rev="v0.1"))
```

## Install CliMA Land Modules (dev)
```julia
using Pkg;
Pkg.add("EmeraldConstants");
Pkg.add("WaterPhysics");
Pkg.add("EmeraldNamespace");
Pkg.add("LeafOptics");
Pkg.add("CanopyRadiativeTransfer");
Pkg.add("Photosynthesis");
Pkg.add("SoilHydraulics");
Pkg.add("PlantHydraulics");
Pkg.add("StomataModels");
Pkg.add("SoilPlantAirContinuum");
```
