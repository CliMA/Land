# CliMA Land

## 关于

CliMA Land 是下一代陆面模式（Land Surface Model, LSM），旨在使用广泛可得的遥感数据以及地面通量观测数据。CliMA Land 是一个高度模块化的平台，用于[...]。

这个项目旨在成为一个社区合作项目，利用世界各地陆面模式研究小组所完成的工作。最终目标是构建一个生物物理[...]  
- 叶片层级的太阳诱导叶绿素荧光（SIF），并通过冠层传播；  
- 各波段的反射率，可由航天器观测到；  
- 土壤和植被的含水量。

CliMA Land 的一个特定关注点是通过测试最新的优化理论以及植物生理学发展来研究水-碳反馈。我们将尽力实现最新的并遵循一些 [c[...]。

整个模型将使用 [Julia](https://docs.julialang.org/en/v1/) 编写（Julia: ["看起来像 Python，感觉像 Lisp，运行像 Fortran"](https://www.youtube.com/watch?v=8h8rQyEpiZA&t=)），其[...]。

![来自太空的荧光](docs/src/assets/world_sif.jpg?raw=true "来自太空的 SIF（太阳诱导荧光）")


## 陆面模式
|||
|:-------------------|:--------------------------------------------|
| **文档构建**       | [![docs build][docs-bld-img]][docs-bld-url] |
| **文档（开发版）** | [![dev][docs-dev-img]][docs-dev-url]        |
| **单元测试**       | [![unit test][st-img]][st-url]              |
| **代码覆盖率**     | [![codecov][codecov-img]][codecov-url]      |
| **Bors**           | [![Bors enabled][bors-img]][bors-url]       |
|||

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


## 示例

### 在单个站点上运行 CliMA Land（v0.1）

有关 CliMA Land v0.1 的示例，请参见 https://github.com/CliMA/Land/tree/v0.1。


## 参考文献

在使用 CliMA Land（v0.1）时，请引用以下文献。

### 一般模型描述
Y. Wang, P. Köhler, L. He, R. K. Braghiere, R. Doughty, J. Wood, C. Frankenberg. 2021.  
Testing stomatal models at the stand level in deciduous angiosperm and evergreen gymnosperm forests using CliMA Land (v0.1).  
Geoscientific Model Development. 14(11): 6741-6763.  
[DOI](https://doi.org/10.5194/gmd-14-6741-2021)  
[PDF](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG/raw/master/publications/wang2021testing.pdf)  
[SI（补充信息）](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG/raw/master/publications/wang2021testing-si.pdf)  
[代码](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG)

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

### 凝集指数（Clumping index）实现
R. K. Braghiere, Y. Wang, R. Doughty, D. Souza, T. Magney, J. Widlowski, M. Longo, A. Bloom, J. Worden, P. Gentine, and C. Frankenberg. 2021.  
Accounting for canopy structure improves hyperspectral radiative transfer and sun-induced chlorophyll fluorescence representations in a new generation Earth System model.  
Remote Sensing of Environment. 261: 112497.  
[DOI](https://doi.org/10.1016/j.rse.2021.112497)  
[PDF](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG/raw/master/publications/braghiere2021accounting.pdf)  
[SI（补充信息）](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG/raw/master/publications/braghiere2021accounting-si.pdf)  
[代码](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG)

```
@article{braghiere2021accounting,
    author = {Braghiere, Renato K and Wang, Yujie and Doughty, Russell and Sousa, Daniel and Magney, Troy and Widlowski, Jean-Luc and Longo, Marcos and Bloom, A Anthony and Worden, John and Gentine, P[...]
    year = {2021},
    title = {Accounting for canopy structure improves hyperspectral radiative transfer and sun-induced chlorophyll fluorescence representations in a new generation Earth System model},
    journal = {Remote Sensing of Environment},
    volume = {261},
    pages = {112497}
}
```

### 冠层复杂性
Y. Wang, C. Frankenberg. 2022.  
On the impact of canopy model complexity on simulated carbon, water, and solar-induced chlorophyll fluorescence fluxes.  
Biogeosciences. 19(1): 29-45.  
[DOI](https://doi.org/10.5194/bg-19-29-2022)  
[PDF](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG/raw/master/publications/wang2022impact.pdf)  
[代码](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG)

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
