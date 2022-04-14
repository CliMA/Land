# CLIMA-Land


## About

Note that CliMA Land is being refactored. During this process, no new feature will be added and only bug fix is allowed. We wish to present the refactored CliMA Land in the near future, please wait for the new version v0.2.

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

1. cd into a folder to start with, for example
```shell
$ mkdir Test-CliMA-Land
$ cd Test-CliMA-Land
```

2. install packages from CliMA Land with version control, we recommend to download our preconfigurated and tested Project.toml
```toml
[deps]
CanopyLayers = "677f5362-5107-42e4-8e81-51d9c4a1f96c"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
GriddingMachine = "f20cf718-bf4d-4727-bc8f-485b1f283ac6"
LazyArtifacts = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
Photosynthesis = "537ec7c9-aaee-45d0-8af0-4b77892958a6"
PkgUtility = "0d262f2c-28e9-492c-8e19-d7a5c4f11611"
PlantHydraulics = "d6acb6ec-f7e4-548d-b108-f8a3d9a6ce13"
SoilPlantAirContinuum = "2d76e174-5bec-4df2-b5ea-844408736dc2"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StomataModels = "5fa394f2-99ce-4e3a-9704-be7b3526889d"
UnPack = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
WaterPhysics = "20dd5ee6-61da-454b-ac5d-c09c2977e03a"

[compat]
CanopyLayers = "0.1.15"
DataFrames = "1.3.2"
GriddingMachine = "0.1.8"
Photosynthesis = "0.2.0"
PkgUtility = "=0.1.13"
PlantHydraulics = "0.2.13"
SoilPlantAirContinuum = "0.1.15"
StomataModels = "0.1.11"
UnPack = "1.0.2"
julia = "1.6"
```

Alternatively, you may download it directly from our FTP
```shell
$ wget ftp://fluo.gps.caltech.edu/data/CliMA/tutorial-v0.1/Project.toml
```

3. initialize the Julia environment
```shell
$ julia --project -e "using Pkg; Pkg.instantiate();"
```

4. download the wrapper functions and file we prepare
```shell
$ wget ftp://fluo.gps.caltech.edu/data/CliMA/tutorial-v0.1/clima-land.jl
$ wget ftp://fluo.gps.caltech.edu/data/CliMA/tutorial-v0.1/era5_2019_117_296_1X.csv
```

5. Make sure you have these files in the folder before heading to next step
   - `Project.toml`
   - `Manifest.toml`
   - `clima-land.jl`
   - `era5_2019_55_329_1X.csv`

6. run the model for a site at 26.5N 115.5E for the year 2019
```shell
$ julia --project
```
```julia
julia> include("clima-land.jl");
julia> params = query_data(26.5, 115.5, 2019);
julia> clima_land!(params...);
```

7. you should get a file named `era5_2019_117_296_1X.simulation.hs.csv` after a few minutes

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
