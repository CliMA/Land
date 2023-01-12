# API
```@meta
CurrentModule = PlantHydraulics
```


## About
PlantHydraulics.jl provides numerical hydraulic models that can be used to simulate the flow and pressure profiles in plants. In the model, hydraulic system is a combination of hydraulic organs, which are Root, Stem, and Leaf. And as a result, the functions to simulate flow and pressure profiles can be used at organ level and system level. For example, for the simplest case, one can simulate the pressure profile within a xylem segment; and for the most complex case, one can simulate the flow and pressure profiles in a tree with multiple roots, trunk, multiple canopy layers (a stem and a leaf per layer).

Although PlantHydraulics.jl is meant to server as a fundamental dependency package of the CliMA Land, it can also be used as a standalone package and thus to use with other vegetation modeling other than CliMA Land. Below, we will show some examples of how to use PlantHydraulics.jl at organ and plant level.


## Organ Level Simulations
Root, Stem, and Leaf organs differ in the following:
- Root has a rhizosphere component before the xylem
- Leaf has an extraxylary component after the xylem
- Stem and Root have a height change to account for gravity
- Stem and Root capacitance is along the flow path
- Leaf capacitance is at the end of the flow path (along with extraxylary component)

For example in a Leaf, to simulate the flow and pressure profiles, what you need to do are
```julia
using EmeraldNamespace
using PlantHydraulics
FT = Float64;

leaf = EmeraldNamespace.Leaf{FT}();
stem = EmeraldNamespace.Stem{FT}();
root = EmeraldNamespace.Root{FT}();
soil = EmeraldNamespace.SoilLayer{FT}();
PlantHydraulics.xylem_flow_profile!(leaf, FT(1));
PlantHydraulics.xylem_flow_profile!(stem, FT(1));
PlantHydraulics.xylem_flow_profile!(root, FT(1));
PlantHydraulics.xylem_pressure_profile!(leaf);
PlantHydraulics.xylem_pressure_profile!(stem);
PlantHydraulics.xylem_pressure_profile!(root, soil);
```

Note the second parameter `FT(1)` is the time step in second to use with non-steady state mode. If the flow mode is steady state, `FT(1)` will not be used; otherwise, water source/sink term will be applied based on the state of the capacitance tissue. As there is a rhizosphere component in the Root, which uses soil moisture retension curve to compute soil water potential, we need to pass soil information to the `xylem_pressure_profile!` function when updating soil pressure profile.


## Plant Level Simulation
While one can always customize the hydraulic system using functions `xylem_flow_profile!` and `xylem_pressure_profile!`, we provide shortcut functions to soil-plant-air continuum where soil layers, plant, and air layers are aligned. For example, a spac of a tree (`MonoMLTreeSPAC`).
```julia
spac = EmeraldNamespace.MonoMLTreeSPAC{FT}();
PlantHydraulics.xylem_flow_profile!(spac, FT(1));
PlantHydraulics.xylem_pressure_profile!(spac);
```

Similarly, one can use `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, and `MonoElementSPAC` (`Mono` means mono species spac, `ML` means multiple root and canopy layers, `Grass` with no trunk and branch, `Palm` with no branch, and `Element` means the spac consists of one soil layer, root, stem, and leaf element each). Currently, the predefined SPAC only supports the four listed above. If you need more customized SPACs, you will need to combine the fundamental functions `xylem_flow_profile!` and `xylem_pressure_profile!` appropriately. Be carefuly to not forget to synchronize the flow rates with stomtal conductance (this is done automactically for predefined SPACs).


## Vulnerability curve
```@docs
relative_hydraulic_conductance(vc::ComplexVC{FT}, p_25::FT) where {FT<:AbstractFloat}
critical_pressure
```

## Pressure volume curve
```@docs
xylem_pressure
capacitance_buffer
```

## Cavitation legacy
```@docs
clear_legacy!
```

## FLow profile
```@docs
flow_in
flow_out
root_pk
xylem_flow_profile!
xylem_flow_profile!(organ::Union{Leaf{FT}, Leaves2D{FT}, Root{FT}, Stem{FT}}, Δt::FT) where {FT<:AbstractFloat}
xylem_flow_profile!(spac::MonoElementSPAC{FT}, Δt::FT) where {FT<:AbstractFloat}
```

## Pressure profile
```@docs
xylem_end_pressure
xylem_pressure_profile!
xylem_pressure_profile!(spac::MonoElementSPAC{FT}; update::Bool = true) where {FT<:AbstractFloat}
∂E∂P
```

## Critical flow
```@docs
critical_flow
critical_flow(hs::LeafHydraulics{FT}, T::FT, ini::FT = FT(0.5); kr::FT = FT(0.001)) where {FT<:AbstractFloat}
critical_flow(spac::MonoElementSPAC{FT}, ini::FT = FT(0.5); kr::FT = FT(0.001)) where {FT<:AbstractFloat}
```

## Tuning factor
```@docs
β_factor
β_factor!
β_factor!(spac::MonoElementSPAC{FT}) where {FT<:AbstractFloat}
β_factor!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat}
```


## Energy budget
```@docs
plant_energy!
plant_energy!(spac::MonoMLGrassSPAC{FT}) where {FT<:AbstractFloat}
plant_energy!(spac::MonoMLGrassSPAC{FT}, δt::FT) where {FT<:AbstractFloat}
```
