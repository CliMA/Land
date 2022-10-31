# Photosynthesis.jl

## Install
```julia
using Pkg;
Pkg.add("Photosynthesis");
```

## About
Photosynthesis models for C3 and C4 photosynthesis. Photosynthesis.jl supports three photosynthesis models and two fluorescence models. The photosynthesis models are
- C3VJP model based on classic C3 model, which is known as FvCB model (Farquhar et al. 1980)
- C4VJP model based on classic C4 model, which is known as Collaz model (Collaz et al. 1992)
- C3Cytochrome model based on a new C3 model developed by Johnson and Berry (2021)

We, however, made some modifications by adding a product limited photosynthetic rate to the C3 models, and a Rubisco limited photosynthetic rate to the C4 model.

Besides the tranditional photosynthesis model, we also included functions to compute fluorescence related parameters, such as fluorescence quantum yield and non-photochemical quenching. he two implemented fluorescence models are
- Van der Tol et al. (2013) fluorescence model to use with C3VJP and C4VJP models
- Johnson and Berry (2021) fluorescence model to use with C3Cytochrome model

We aim to make Photosynthesis.jl a standalone package rather than just part of the CliMA Land model. Thus, in the documentations below, we will present examples of how to use Photosynthesis.jl at the leaf level. Same logic applies to canopy scale simulations.


## Model Selection
Starting from v0.3, photosynthesis and fluorescence model selection is done by setting up the fields of a leaf. There are three types of leaf in Photosynthesis (all the structures are defined in ClimaCache.jl and shared among all CliMA Land submodules), and they are
- `Leaf` for a single leaf to use in leaf level research
- `Leaves1D` for a vector of leaves to use in big leaf models
- `Leaves2D` for a matrix of sunlit fractions and a shaded fraction to use along with canopy with leaf angular distribution

For all of the three leaf structs, there are two fields named
- `PSM` for photosynthesis model
- `PRC` for photosynthesis reaction center

A `C3VJPModel` type `PSM` along with a `VJPReactionCenter` type `PRC` defines the C3VJP model; a `C4VJPModel` type `PSM` along with a `VJPReactionCenter` type `PRC` defines the C4VJP model; and a `C3CytochromeModel` type `PSM` along with a `CytochromeReactionCenter` type `PRC` defines the C3Cytochrome model. For instance, `leaf_c3`, `leaf_c4`, and `leaf_cy` each defines a model to use the three predefined photosynthesis models:
```julia
using ClimaCache;
FT = Float64;

leaf_c3 = ClimaCache.Leaf{FT}();
leaf_c4 = ClimaCache.Leaf{FT}(PSM = ClimaCache.C4VJPModel{FT}());
leaf_cy = ClimaCache.Leaf{FT}(PSM = ClimaCache.C3CytochromeModel{FT}(), PRC = ClimaCache.CytochromeReactionCenter{FT}());

# users can define the same fields for Leaves1D and Leaves2D to custoimize photosynthesis model
leaf_d3 = ClimaCache.Leaves1D{FT}();
leaf_d4 = ClimaCache.Leaves1D{FT}(PSM = ClimaCache.C4VJPModel{FT}());
leaf_dy = ClimaCache.Leaves1D{FT}(PSM = ClimaCache.C3CytochromeModel{FT}(), PRC = ClimaCache.CytochromeReactionCenter{FT}());
leaf_e3 = ClimaCache.Leaves2D{FT}();
leaf_e4 = ClimaCache.Leaves2D{FT}(PSM = ClimaCache.C4VJPModel{FT}());
leaf_ey = ClimaCache.Leaves2D{FT}(PSM = ClimaCache.C3CytochromeModel{FT}(), PRC = ClimaCache.CytochromeReactionCenter{FT}());
```


## Photosynthesis Model Procedure
For all three photosynthesis+fluorescence combo models, photosynthetic rates are computed using the following procedure:
- Update temperature dependent variables using `photosystem_temperature_dependence!`
- Calculate electron transport rate using `photosystem_electron_transport!`
- Calculate RubisCO limited rate using `rubisco_limited_rate!`
- Calculate light limited rate using `light_limited_rate!`
- Calculate product limited rate using `product_limited_rate!`
- Calculate gross and net rates using `colimit_photosynthesis!`
- Update fluorescence related variables using `photosystem_coefficients!`

Yet, for convenience, all the listed steps are combined in one function `leaf_photosynthesis!` (see page API for all supported methods). At leaf level, one can simply call the function, for example
```julia
using Photosynthesis

air    = ClimaCache.AirLayer{FT}();
g_mode = ClimaCache.GCO₂Mode();
p_mode = ClimaCache.PCO₂Mode();
Photosynthesis.leaf_photosynthesis!(leaf_c3, air, g_mode);
Photosynthesis.leaf_photosynthesis!(leaf_cy, air, p_mode);
```

Note here that (1) we need parameters from surrounding air to compute photosynthetic rates, such as oxygen concentration, (2) `GCO₂Mode` mode is used when leaf total diffusive conductance for CO₂ is known, and `PCO₂Mode` is used when leaf internal CO₂ partial pressure is known, and (3) leaf and air conditions must be synced before calling the `leaf_photosynthesis!` function. For example, to construct an A-Ci curve, one will need to change the field `_p_CO₂_i` to proceed:
```julia
leaf_c3.t = 300;
leaf_c3.ppar = 1000;
for _p in 5:5:100
    leaf_c3._p_CO₂_i = _p;
    Photosynthesis.leaf_photosynthesis!(leaf_c3, air, p_mode);
    @info "Photosynthetic rate at" leaf_c3._p_CO₂_i leaf_c3.a_net;
end
```

In the example above, we first defined leaf temperature is 300 `K` (note that we use Kelvin here, not Celcius), then prescribed a PPAR of 1000 `μmol m⁻² s⁻¹`. PPAR is the photosynthetically active radiation (PAR) that goes to photosystems; PPAR is different from APAR (PAR absorbed by leaf, some APAR does not go to photosystems). Ssee our [technical note](
https://doi.org/10.5194/bg-2022-172) for the difference between PAR, APAR, and PPAR. Similarly, we can try out how Anet responds to leaf temperature and PPAR by modifying the example above.
