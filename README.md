# CLIMA-Land
This project is supposed to be a community effort, leveraging all the work that has been done in Land Surface Modeling from various groups around the world. The ultimate goal here is to build a Soil-Plant-Atmosphere continuum (SPAC) bio-physical model that represents the state of the art and can be coupled to the Climate-Machine, i.e. Caltech's CLIMA initiative. The model will include multi-layer soil and canopy properties and models "observables" that can be used as constraint, e.g. Solar Induced Chlorophyll Fluorescence (SIF) on the leaf-level propagated through the canopy, reflectance in various bands as measured from space, soil and vegetation moisture content. A specific focus will be on water-carbon feedbacks by testing recent developments in stomatal optimization theories as well as plant hydraulics. We will try to adhere to some [coding structure](https://github.com/gbonan/CLM-ml_v0) developed by Gordan Bonan but implement parts from other programs, such as [SCOPE](https://github.com/Christiaanvandertol/SCOPE).

The entire model will be written in [Julia](https://docs.julialang.org/en/v1/) (Julia: ["Looks like Python, feels like Lisp, runs like Fortran"](https://www.youtube.com/watch?v=8h8rQyEpiZA&t=)), which should make the barrier of entry lower for incoming students, PostDocs, etc). If you want to contribute, please contact us (cfranken@caltech.edu).

![Fluorescence from Space](pics/world_sif.jpg?raw=true "SIF from Space")

|||
|---------------------:|:----------------------------------------------|
| **Docs Build**       | [![docs build][docs-bld-img]][docs-bld-url]   |
| **Documentation**    | [![dev][docs-dev-img]][docs-dev-url]          |
| **Azure Build**      | [![azure][azure-img]][azure-url]              |
| **Code Coverage**    | [![codecov][codecov-img]][codecov-url]        |
| **Bors**             | [![Bors enabled][bors-img]][bors-url]         |

[docs-bld-img]: https://github.com/climate-machine/Land/workflows/Documentation/badge.svg
[docs-bld-url]: https://github.com/climate-machine/Land/actions?query=workflow%3ADocumentation

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://climate-machine.github.io/Land/dev/

[azure-img]: https://dev.azure.com/climate-machine/Land/_apis/build/status/climate-machine.Land?branchName=master
[azure-url]: https://dev.azure.com/climate-machine/Land/_build/latest?definitionId=1&branchName=master

[codecov-img]: https://codecov.io/gh/climate-machine/Land/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/climate-machine/Land

[bors-img]: https://bors.tech/images/badge_small.svg
[bors-url]: https://app.bors.tech/repositories/12357
