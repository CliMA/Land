# CLIMA-Land


## About

Note that CliMA Land is being refactored to v0.2. As a result, no new feature will be added and only bug fixes are allowed. We wish to present the refactored CliMA Land in the near future, please wait for the new version v0.2.

This project is supposed to be a community effort, leveraging all the work that has been done in Land Surface Modeling from various groups around the world. The ultimate goal here is to build a Soil-Plant-Atmosphere continuum (SPAC) bio-physical model that represents the state of the art and can be coupled to the Climate-Machine, i.e. Caltech's CLIMA initiative. The model will include multi-layer soil and canopy properties and models "observables" that can be used as constraint, e.g. Solar Induced Chlorophyll Fluorescence (SIF) on the leaf-level propagated through the canopy, reflectance in various bands as measured from space, soil and vegetation moisture content. A specific focus will be on water-carbon feedbacks by testing recent developments in stomatal optimization theories as well as plant hydraulics. We will try to adhere to some [coding structure](https://github.com/gbonan/CLM-ml_v0) developed by Gordan Bonan but implement parts from other programs, such as [SCOPE](https://github.com/Christiaanvandertol/SCOPE).

The entire model will be written in [Julia](https://docs.julialang.org/en/v1/) (Julia: ["Looks like Python, feels like Lisp, runs like Fortran"](https://www.youtube.com/watch?v=8h8rQyEpiZA&t=)), which should make the barrier of entry lower for incoming students, PostDocs, etc). If you want to contribute, please contact us (cfranken@caltech.edu).

![Fluorescence from Space](pics/world_sif.jpg?raw=true "SIF from Space")


## Land model
|||
|:-------------------|:--------------------------------------------|
| **Docs Build**     | [![docs build][docs-bld-img]][docs-bld-url] |
| **Documentation**  | [![dev][docs-dev-img]][docs-dev-url]        |
| **Unit Test**      | [![unit test][st-img]][st-url]              |
| **Code Coverage**  | [![codecov][codecov-img]][codecov-url]      |
| **Bors**           | [![Bors enabled][bors-img]][bors-url]       |

[docs-bld-img]: https://github.com/CliMA/Land/workflows/Documentation/badge.svg
[docs-bld-url]: https://github.com/CliMA/Land/actions?query=workflow%3ADocumentation

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://CliMA.github.io/Land/dev/

[st-img]: https://github.com/CliMA/Land/workflows/JuliaStable/badge.svg?branch=main
[st-url]: https://github.com/CliMA/Land/actions?query=workflow%3AJuliaStable

[codecov-img]: https://codecov.io/gh/CliMA/Land/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/Land

[bors-img]: https://bors.tech/images/badge_small.svg
[bors-url]: https://app.bors.tech/repositories/24777


## Examples

### Run CliMA Land for a single site (v0.1)

1. Download all the files in folder examples into a folder,
```shell
$ mkdir Test-CliMA-Land
$ cd Test-CliMA-Land
$ wget https://github.com/CliMA/Land/raw/v0.1/examples/debug.jld2
$ wget https://github.com/CliMA/Land/raw/v0.1/examples/debug.nc
$ wget https://github.com/CliMA/Land/raw/v0.1/examples/example.jl
$ wget https://github.com/CliMA/Land/raw/v0.1/examples/Manifest.toml
$ wget https://github.com/CliMA/Land/raw/v0.1/examples/Project.toml
```

2. initialize the Julia environment
```shell
$ julia --project -e "using Pkg; Pkg.instantiate();"
```

3. Make sure you have these files in the folder before heading to next step
   - `Project.toml`
   - `Manifest.toml`
   - `example.jl`
   - `debug.jld2`
   - `debug.nc`

4. run the model for the site
```shell
$ julia --project example.jl
```

5. you should get a file named `debug.output.nc` after a few minutes

You will need to edit the functions we provided if you need to change more parameters, or output more results. Feel free to contact us through email or Github Issues (if this tutorial does not work).


## References

Please cite the following when you use the CliMA Land (v0.1)

### General model description
Y. Wang, P. Köhler, L. He, R. K. Braghiere, R. Doughty, J. Wood, C. Frankenberg. 2021.
Testing stomatal models at the stand level in deciduous angiosperm and evergreen gymnosperm forests using CliMA Land (v0.1).
Geoscientific Model Development. 14(11): 6741-6763.
[DOI](https://doi.org/10.5194/gmd-14-6741-2021)
[PDF](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG/raw/master/publications/wang2021testing.pdf)
[SI](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG/raw/master/publications/wang2021testing-si.pdf)
[CODE](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG)

```
@article{wang2021testing,
    author = {Wang, Y. and K{\"o}hler, P. and He, L. and Doughty, R. and Braghiere, R. K. and Wood, J. D. and Frankenberg, C.},
    year = {2021},
    title = {Testing stomatal models at the stand level in deciduous angiosperm and evergreen gymnosperm forests using CliMA Land (v0.1)},
    journal = {Geoscientific Model Development},
    volume = {14},
    number = {11},
    pages = {6741--6763}
}
```

### Clumping index implementation
R. K. Braghiere, Y. Wang, R. Doughty, D. Souza, T. Magney, J. Widlowski, M. Longo, A. Bloom, J. Worden, P. Gentine, and C. Frankenberg. 2021.
Accounting for canopy structure improves hyperspectral radiative transfer and sun-induced chlorophyll fluorescence representations in a new generation Earth System model.
Remote Sensing of Environment. 261: 112497.
[DOI](https://doi.org/10.1016/j.rse.2021.112497)
[PDF](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG/raw/master/publications/braghiere2021accounting.pdf)
[SI](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG/raw/master/publications/braghiere2021accounting-si.pdf)
[CODE](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG)

```
@article{braghiere2021accounting,
    author = {Braghiere, Renato K and Wang, Yujie and Doughty, Russell and Sousa, Daniel and Magney, Troy and Widlowski, Jean-Luc and Longo, Marcos and Bloom, A Anthony and Worden, John and Gentine, Pierre and Frankenberg, Christian},
    year = {2021},
    title = {Accounting for canopy structure improves hyperspectral radiative transfer and sun-induced chlorophyll fluorescence representations in a new generation Earth System model},
    journal = {Remote Sensing of Environment},
    volume = {261},
    pages = {112497}
}
```

### Canopy complexity
Y. Wang, C. Frankenberg. 2022.
On the impact of canopy model complexity on simulated carbon, water, and solar-induced chlorophyll fluorescence fluxes.
Biogeosciences. 19(1): 29-45.
[DOI](https://doi.org/10.5194/bg-19-29-2022)
[PDF](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG/raw/master/publications/wang2022impact.pdf)
[CODE](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG)

```
 @article{wang2022impact,
 	author = {Wang, Y. and Frankenberg, C.},
 	year = {2022},
 	title = {On the impact of canopy model complexity on simulated carbon, water, and solar-induced chlorophyll fluorescence fluxes},
 	journal = {Biogeosciences},
 	volume = {19},
 	number = {1},
 	pages = {29--45}
 }
 ```
