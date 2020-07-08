# API
```@meta
CurrentModule = PlantHydraulics
```




## Plant Hydraulic System
The PlantHydraulics module provides two levels of hydraulics system: organ-level and plant-level. The organ-level hydraulic systems include Leaf, Root, and Stem (trunk and branch). The plant-level hydraulic system is can be any combination of the three organs (custimized definition may apply).




## Leaf, Root, and Stem HS
Plant hydraulics is segmented to three organ-level systems/structs ([`LeafHydraulics`](@ref), [`RootHydraulics`](@ref), and [`StemHydraulics`](@ref)) subject to an Abstract type ([`AbstractHydraulicSystem`](@ref)). The major differences among the three structs are
- [`LeafHydraulics`](@ref) has an extra-xylary component
- [`RootHydraulics`](@ref) has a rhizosphere component
- [`RootHydraulics`](@ref) and [`StemHydraulics`](@ref) have a gravity component

See the documentation for each struct for more details:
```@docs
AbstractHydraulicSystem
LeafHydraulics
RootHydraulics
StemHydraulics
```

To initialize a hydraulics system, one needs to provide the floating type, for example:
```julia
using PlantHydraulics

FT = Float32;
hs_leaf = LeafHydraulics{FT}();
hs_root = RootHydraulics{FT}();
hs_stem = StemHydraulics{FT}();
```




## Whole-plant HS
Plants differ in their structures, for example, some plants have a canopy far above the ground elevated by a trunk, some plants have a structured canopy supported by branch systems, and some plant has no trunk at all. To represent the structural differences, several types of plant hydraulics systems are pre-defined, and they are [`GrassLikeHS`](@ref), [`PalmLikeHS`](@ref), and [`TreeLikeHS`](@ref) structs subject to a [`AbstractPlantHS`](@ref) type, where `HS` stands for hydraulic system. The major difference between the `HS`s are
- [`GrassLikeHS`](@ref) has only mutiple root and canopy layers, no trunk or branch
- [`PalmLikeHS`](@ref) has multiple root layers, a trunk, and multiple canopy layers, no branch system
- [`TreeLikeHS`](@ref) has multiple root layers, a trunk, and multiple branch+canopy layers, and each branch corresponds to a canopy layer

See the documentation for each struct for more details:
```@docs
AbstractPlantHS
GrassLikeHS
PalmLikeHS
TreeLikeHS
```

To ease the initialization of a plant hydraulics system, a few customized functions are provided for quick initialization. More importantly, modifications to each field in the struct are always allowed. The quick functions are [`create_grass_like_hs`](@ref), [`create_palm_like_hs`](@ref), and [`create_tree_like_hs`](@ref):
```@docs
create_grass_like_hs
create_palm_like_hs
create_tree_like_hs
```

What these functions do are to determine how many root layers and branch/canopy layers to add based on the tree information and environmental settings. To determine number of root layers, rooting depth and the soil layer information are required. The `z_root` is the maximal root depth in negative number, and `soil_bounds` is the boundaries of soil layers staring from 0. For example, for a `soil_bounds` of [0.0, -1.0, -2.0, -3.0, -4.0], a `z_root` of -1 gives 1 root layer, and a `z_root` of -1.5 or -2.0 gives 2 root layers. The `z_trunk`, `z_canopy`, and `air_bounds` determine how many canopy layers to add. For example, for a `air_bounds` of [0.0, 1.0, 2.0, 3.0, 4.0, 5.0 ... 20.0, 21.0, 22.0], a `z_trunk` of 5.0 `z_canopy` of 7.0 give 2 canopy layers, and a `z_trunk` of 5.5 `z_canopy` of 7.0 give 2 canopy layers. Also, the `root_index_in_soil` and `canopy_index_in_air` indicate which soil or air layer the root or canopy layer corresponds with, respectively. For instance, a index of 7 means that the canopy layer should use the 7th air layer.

To initialize a whole-plant hydraulic system, checkout the example below:
```julia
using PlantHydraulics

FT = Float32;
grass = create_grass_like_hs(FT(-2.1), FT(0.5), FT(8), FT[0,-1,-2,-3], collect(FT,0:1:20));
palm  =  create_palm_like_hs(FT(-2.1), FT(0.5), FT(8), FT[0,-1,-2,-3], collect(FT,0:1:20));
tree  =  create_tree_like_hs(FT(-2.1), FT(0.5), FT(8), FT[0,-1,-2,-3], collect(FT,0:1:20));
```




## Hydraulic conductance
Plants transport water through xylem conduits (vessels in most angiosperms, trachieds in most gymnosperms). With the ascent of sap along the hydraulic system, water pressure in the conduits is typically negative. The negative xylem water pressure tends to pull air from surrounding tisses or the atmosphere into the xylem conduits, resulting in xylem cavitation. The air bubbles in cavitated conduits block water flow, and thus results in decline of water transport capability (measured by xylem hydraulic conductance).

Typically, the correlation between xylem water pressure ($P \leq 0$) and hydraulic conductance ($k$) is expressed by a Weibull function for [`WeibullSingle`](@ref) type correlation:
```math
k = k_\text{max} \cdot \exp \left( -\left( \dfrac{-P}{B} \right)^C \right)
```
where $k_\text{max}$ is the maximal hydraulic conductance, and $B$ and $C$ are the Weibull parameters. This correlation is also known as vulnerability curve (VC) to drought stress. Sometimes, plants exhibit a segmented VC, for example, the fibers may transport water as well and are much more resistant to drought than vessels. Thus, a dual Weibull function is presented for [`WeibullDual`](@ref) type correlation ($P \leq 0$):
```math
k = k_\text{max} \cdot \left\{ f_1 \cdot \exp \left[ -\left( \dfrac{-P}{B_1} \right)^{C_1} \right] +
                         (1 - f_1) \cdot \exp \left[ -\left( \dfrac{-P}{B_2} \right)^{C_2} \right] \right\}
```

The VC formulations are abstractized as
```@docs
AbstractVulnerability
WeibullDual
WeibullSingle
```

The function to call is
```@docs
xylem_k_ratio
```
Note it here that `xylem_k_ratio(vc, p)` calculate the k without making temperature corrections, but the `xylem_k_ratio(vc, p_25, vis)` makes correction over the viscosity (the higher the viscosity, the lower the k). Also, `p_25` means that the pressure has been corrected to 298.15 K for surface tension (the higher the surface tension, the more resistant the xylem).

Meanwhile, there is a function to call to calculate the critical pressure, beyond which leaf will decicate. The critical pressure is calculated as the pressure at which $k$ is 0.001 of $k_\text{max}$ for [`WeibullSingle`](@ref) (for [`WeibullDual`](@ref), each segment need to reach 0.001). The functions is
```@docs
xylem_p_crit
```

Examples:
```julia
using PlantHydraulics

FT = Float32;
vc_1 = WeibullSingle{FT}();
vc_2 = WeibullDual{FT}();

k_1 = xylem_k_ratio(vc_1, -1.0);
k_2 = xylem_k_ratio(vc_2, -1.0);
k_3 = xylem_k_ratio(vc_1, -1.0, 1.2);
k_4 = xylem_k_ratio(vc_2, -1.0, 1.2);
```




## Pressure and Flow
The PlantHydraulics module is designed to run numerically for the following reasons:
- Weibull function is cannot be integrated
- The VC is segmented, i.e., if $P > 0$, $k = k_\text{max}$ (implemented in [`xylem_k_ratio`](@ref))
- Once xylem cavitation occurs, it cannot be easily recovered unless $P > 0$, and thus there is a drought legacy effect. This is why there are a few fields in the [`LeafHydraulics`](@ref), [`RootHydraulics`](@ref), and [`StemHydraulics`](@ref) structs to store the drought history information.
- Temperature may change along the flow path. The `f_st` and `f_vis` in the structs help deal with these effects.

Function [`xylem_p_from_flow`](@ref) calculates the xylem end pressure for an organ-level hysraulic system. As mentioned above, the [`RootHydraulics`](@ref) and [`StemHydraulics`](@ref) has a gravity component, and the [`RootHydraulics`](@ref) has a rhizosphere component. Also be aware that [`xylem_p_from_flow`](@ref) accounts for temperature effects on surface tension and viscosity.
```@docs
xylem_p_from_flow
```

Noe that function [`xylem_p_from_flow`](@ref) does not update the pressure profiles or history in the xylem. To update these profiles, use [`hydraulic_p_profile!`](@ref):
```@docs
hydraulic_p_profile!
```

Examples:
```julia
using PlantHydraulics

FT = Float32;
leaf = LeafHydraulics{FT}();
p = xylem_p_from_flow(leaf, FT(0.01));
@show leaf.p_element;
hydraulic_p_profile!(leaf, FT(0.01));
@show leaf.p_element;
```




## Root Hydraulics
Function [`xylem_p_from_flow`](@ref) works for the case of only 1 root layer if one needs the plant base xylem water pressure. However, when there are multiple root layers, [`xylem_p_from_flow`](@ref) does not apply. In this case, iterations are required to calculate the xylem end pressure for each root layers, and then make sure all root layers have the same xylem end pressure. A few functions are provided to realize this.

Function [`root_q_from_pressure`](@ref) uses Root Solving method to calculate the flow rate through the [`RootHydraulics`](@ref) struct that yields the given xylem end pressure. The `ini` in the function is optional. However, using the flow rate from last instant when pressure does not differ much will speed up the calculation.
```@docs
root_q_from_pressure
```

Computing the flow rates in all the root layers and the plant base xylem pressure from a given flow needs more iterations based on [`root_q_from_pressure`](@ref). The algorithm is to calculate the flow rates in all the root layers for a given plant base xylem pressure, and then use the Root Solving methods to calculate pressure that yields a total flow rate equals to the given value. The function [`q_diff`](@ref) calculates the difference of total flow rate for a given `pressure` and the `target` flow rate. The function [`root_qs_p_from_q`](@ref) calculates the flows rates in all the root layers and the plant base xylem pressure. Note that `ini` in the function is also optional, but a good guess will speed up the calculations. Be aware that function [`root_qs_p_from_q`](@ref) does not update the pressure profiles in the plant struct.
```@docs
q_diff
root_qs_p_from_q
```

Example
```julia
using PlantHydraulics

FT = Float32;
palm  =  create_palm_like_hs(FT(-2.1), FT(0.5), FT(8), FT[0,-1,-2,-3], collect(FT,0:1:20));
qs,p = root_qs_p_from_q(palm.roots, FT(1));
```




## Leaf Hydraulics
The stomatal models often require plant hydraulics either as a correction factor (in empirical stomatal models) or as the risk term (in optimal stomatal models). To facilitate the calculations, a few specific functions are provided.

Function [`leaf_xylem_risk`](@ref) returns the risk in xylem hydraulic function based on the most downstream end of the xylem. The risk of plant hydraulic system is not only on current system, but also potential new growth (plants don't want to risk new growth either). Thus, function [`leaf_xylem_risk`](@ref) evaluates the risk from the xylem pressure calculated from current system (with drought history), and then compute the risk from the pressure (the severer the srought history, the higher the risk):
```@docs
leaf_xylem_risk
```

Note that function [`leaf_xylem_risk`](@ref) can work on its own without having other organ-level components. For example, by changing the `p_ups` of a [`LeafHydraulics`](@ref), one can simulate the case of drought without caring about other hydraulic systems. Same for function [`leaf_e_crit`](@ref) below. However, these functions are only useful for sensitivity analysis or when `p_ups` in the [`LeafHydraulics`](@ref) is accurate.

Examples
```julia
using PlantHydraulics

FT = Float32;
leaf = LeafHydraulics{FT}();
risk = leaf_xylem_risk(leaf, FT(0.01));
@show risk;
leaf.p_ups = FT(-1.0);
risk = leaf_xylem_risk(leaf, FT(0.01));
@show risk;
```

Function [`leaf_e_crit`](@ref) calculates critical leaf transpiration rate, beyond which leaf will desicate. Function [`leaf_e_crit`](@ref) accounts for drought legacy effect by design, and the more severe the drought history, the lower the `e_crit`. Again, `ini` in the function is also optional, but a good guess will speed up the calculations.
```@docs
leaf_e_crit
```

Examples
```julia
using PlantHydraulics

FT = Float32;
leaf = LeafHydraulics{FT}();
risk = leaf_e_crit(leaf);
@show risk;
leaf.p_ups = FT(-1.0);
risk = leaf_e_crit(leaf);
@show risk;
```

## Whow-plant Hydraulics
Though [`leaf_xylem_risk`](@ref) and [`leaf_e_crit`](@ref) can work on their own, the functions only evaluate the risks on leaf level. The more realistic case is that when leaf transpiration rate increases, `p_ups` in the [`LeafHydraulics`](@ref) gets more negative. Thus, the [`leaf_xylem_risk`](@ref) and [`leaf_e_crit`](@ref) tends to underestimate the risk and overestimate the critical flow rate. To overcome this problem, whole-plant level plant hydraulics are provided.

Function [`tree_p_from_flow`](@ref) calculates the leaf xylem end pressure for a whole-plant struct using these steps:
- calculate the plant base pressure from a given total flow rate
- calculate the trunk end pressure (if present)
- calculate the branch end pressure (if present)
- calculate the leaf end pressure (if present)
```@docs
tree_p_from_flow
```

Accordingly, there is a function [`tree_e_crit`](@ref) to calculate the critical flow rate for the whole plant. Be aware that function [`tree_p_from_flow`](@ref) and [`tree_e_crit`](@ref) only applies to the case of only one canopy layer (or big-leaf model). As to the case of multiple canopy layer, more functions are pending.
```@docs
tree_e_crit
```
