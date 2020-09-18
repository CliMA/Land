# Photosynthesis.jl

<!-- Links and shortcuts -->
[ps-url]: https://github.com/Yujie-W/Photosynthesis.jl
[ps-api]: https://yujie-w.github.io/Photosynthesis.jl/stable/API/

[dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[dev-url]: https://Yujie-W.github.io/Photosynthesis.jl/dev/

[rel-img]: https://img.shields.io/badge/docs-stable-blue.svg
[rel-url]: https://Yujie-W.github.io/Photosynthesis.jl/stable/

[st-img]: https://github.com/Yujie-W/Photosynthesis.jl/workflows/JuliaStable/badge.svg?branch=master
[st-url]: https://github.com/Yujie-W/Photosynthesis.jl/actions?query=branch%3A"master"++workflow%3A"JuliaStable"

[bm-img]: https://github.com/Yujie-W/Photosynthesis.jl/workflows/Benchmarks/badge.svg?branch=master
[bm-url]: https://github.com/Yujie-W/Photosynthesis.jl/actions?query=branch%3A"master"++workflow%3A"Benchmarks"

[min-img]: https://github.com/Yujie-W/Photosynthesis.jl/workflows/Julia-1.3/badge.svg?branch=master
[min-url]: https://github.com/Yujie-W/Photosynthesis.jl/actions?query=branch%3A"master"++workflow%3A"Julia-1.3"

[cov-img]: https://codecov.io/gh/Yujie-W/Photosynthesis.jl/branch/master/graph/badge.svg
[cov-url]: https://codecov.io/gh/Yujie-W/Photosynthesis.jl

## About

Photosynthesis models for C3 and C4 photosynthesis.

| Documentation                                   | CI Status             | Benchmarks            | Compatibility           | Code Coverage           |
|:-----------------------------------------------:|:---------------------:|:---------------------:|:------------------------|:------------------------|
| [![][dev-img]][dev-url] [![][rel-img]][rel-url] | [![][st-img]][st-url] | [![][bm-img]][bm-url] | [![][min-img]][min-url] | [![][cov-img]][cov-url] |




## Dependencies

| Dependency          | Version  | Requirements |
|:--------------------|:---------|:-------------|
| BenchmarkTools      | 0.5.0 +  | Julia 1.0 +  |
| CLIMAParameters     | 0.1.6 +  | Julia 1.3 +  |
| DocStringExtensions | 0.8.3 +  | Julia 0.7 +  |
| Parameters          | 0.12.1 + | Julia 1.0 +  |
| WaterPhysics        | 0.1.0 +  | Julia 1.3 +  |




## Installation
```julia
julia> using Pkg;
julia> Pkg.add("Photosynthesis");
```




## API
See [`API`][ps-api] for more detailed information about how to use [`Photosynthesis.jl`][ps-url].
