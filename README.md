# Photosynthesis.jl

<!-- Links and shortcuts -->
[ps-url]: https://github.com/Yujie-W/Photosynthesis.jl
[ps-api]: https://yujie-w.github.io/Photosynthesis.jl/stable/API/

[dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[dev-url]: https://Yujie-W.github.io/Photosynthesis.jl/dev/

[rel-img]: https://img.shields.io/badge/docs-stable-blue.svg
[rel-url]: https://Yujie-W.github.io/Photosynthesis.jl/stable/

[st-img]: https://github.com/Yujie-W/Photosynthesis.jl/workflows/JuliaStable/badge.svg?branch=main
[st-url]: https://github.com/Yujie-W/Photosynthesis.jl/actions?query=branch%3A"main"++workflow%3A"JuliaStable"

[min-img]: https://github.com/Yujie-W/Photosynthesis.jl/workflows/Julia-1.6/badge.svg?branch=main
[min-url]: https://github.com/Yujie-W/Photosynthesis.jl/actions?query=branch%3A"main"++workflow%3A"Julia-1.6"

[cov-img]: https://codecov.io/gh/Yujie-W/Photosynthesis.jl/branch/main/graph/badge.svg
[cov-url]: https://codecov.io/gh/Yujie-W/Photosynthesis.jl


## Note

`Photosynthesis.jl` has been refactored based on the ClimaCache.jl package. Functions has been renamed to be more specific, but it is getting more user-friendly. As `Photosynthesis.jl` is a module of
    CliMA Land, please cite [CliMA Land](https://github.com/CliMA/Land) when you use `Photosynthesis.jl`.


## About

Photosynthesis models for modeling C3 and C4 photosynthesis.

| Documentation                                   | CI Status             | Compatibility           | Code Coverage           |
|:-----------------------------------------------:|:---------------------:|:------------------------|:------------------------|
| [![][dev-img]][dev-url] [![][rel-img]][rel-url] | [![][st-img]][st-url] | [![][min-img]][min-url] | [![][cov-img]][cov-url] |


## Installation
```julia
julia> using Pkg;
julia> Pkg.add("Photosynthesis");
```


## API
See [`API`][ps-api] for more detailed information about how to use [`Photosynthesis.jl`][ps-url].
