# WaterPhysics.jl

<!-- Links and shortcuts -->
[wp-url]: https://github.com/Yujie-W/WaterPhysics.jl
[wp-api]: https://yujie-w.github.io/WaterPhysics.jl/stable/API/
[cp-url]: https://github.com/CliMA/CLIMAParameters.jl

[dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[dev-url]: https://Yujie-W.github.io/WaterPhysics.jl/dev/

[rel-img]: https://img.shields.io/badge/docs-stable-blue.svg
[rel-url]: https://Yujie-W.github.io/WaterPhysics.jl/stable/

[st-img]: https://github.com/Yujie-W/WaterPhysics.jl/workflows/JuliaStable/badge.svg?branch=master
[st-url]: https://github.com/Yujie-W/WaterPhysics.jl/actions?query=branch%3A"master"++workflow%3A"JuliaStable"

[bm-img]: https://github.com/Yujie-W/WaterPhysics.jl/workflows/Benchmarks/badge.svg?branch=master
[bm-url]: https://github.com/Yujie-W/WaterPhysics.jl/actions?query=branch%3A"master"++workflow%3A"Benchmarks"

[v13-img]: https://github.com/Yujie-W/WaterPhysics.jl/workflows/Julia-1.3/badge.svg?branch=master
[v13-url]: https://github.com/Yujie-W/WaterPhysics.jl/actions?query=branch%3A"master"++workflow%3A"Julia-1.3"

[v14-img]: https://github.com/Yujie-W/WaterPhysics.jl/workflows/Julia-1.4/badge.svg?branch=master
[v14-url]: https://github.com/Yujie-W/WaterPhysics.jl/actions?query=branch%3A"master"++workflow%3A"Julia-1.4"


## About

[`WaterPhysics.jl`][wp-url] includes a collection of temperature dependencies of physical properties of water. Due to its dependency on [`CLIMAParameters.jl`][cp-url], [`WaterPhysics.jl`][wp-url] only supports Julia 1.3 and above.

| Documentation                                   | CI Status             | Benchmarks            | Compatibility                                   |
|:------------------------------------------------|:----------------------|:----------------------|:------------------------------------------------|
| [![][dev-img]][dev-url] [![][rel-img]][rel-url] | [![][st-img]][st-url] | [![][bm-img]][bm-url] | [![][v14-img]][v14-url] [![][v13-img]][v13-url] |




## Dependencies

| Dependency      | Version | Requirements |
|:----------------|:--------|:-------------|
| BenchmarkTools  | 0.5.0 + | Julia 1.0 +  |
| CLIMAParameters | 0.1.6 + | Julia 1.3 +  |




## Installation
```julia
julia> using Pkg;
julia> Pkg.add("WaterPhysics");
```




## API
See [`API`](wp-api) for more detailed information about how to use [`WaterPhysics.jl`][wp-url].
