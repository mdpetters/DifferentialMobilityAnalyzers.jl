# DifferentialMobilityAnalyzers.jl

*A Julia package for working with differential mobility analyzers.*

| **Documentation**                                                               | **Build Status**                                                                                | **Citations** |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] | [![travis badge][travis_badge]][travis_url] [![codecov badge][codecov_badge]][codecov_url] | [![DOI](https://img.shields.io/badge/DOI-10.1080/02786826.2018.1530724-blue?label=DOI)](https://doi.org/10.1080/02786826.2018.1530724) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2652893.svg)](https://doi.org/10.5281/zenodo.2652893)  |

# Installation

The package can be installed from the Julia package prompt with

```julia
julia> ]add  https://github.com/mdpetters/DifferentialMobilityAnalyzers.jl.git
```

The closing square bracket switches to the package manager interface and the ```add``` command installs the package and any missing dependencies. To return to the Julia REPL hit the ```delete``` key.

To load the package run

```julia
julia> using DifferentialMobilityAnalyzers
```


## Documentation

- [**STABLE**][docs-stable-url] &mdash; **documentation of the most recently tagged version.**
- [**DEVEL**][docs-dev-url] &mdash; *documentation of the in-development version.*


## Project Status
The current version of the package is being developed for, Julia `1.4` and above on Linux. It very likely works on macOS and Windows.

The original version v1.0.0 was developed for Julia v0.6. A virtual machine with the original code is archived on zendo. Support for Julia v0.6 was dropped in version 2. 

## Contribute
Contributions including notebooks for classroom instruction, homework assignments, interesting DMA configurations, new inversion schemes, and improved or new functionalities of the language are welcome.

## Citations
This work was supported by the United States Department of Energy, Office of Science, Biological and Environment Research, Grant number DE-SC0018265.

Petters, M.D. (2018) <i> A language to simplify computation of differential mobility analyzer response functions </i> Aerosol Science & Technology, 52 (12), 1437-1451, https://doi.org/10.1080/02786826.2018.1530724.

Petters, M.D. (2019, April 27) <i> Virtual Machine containing Software for "A language to simplify computation of differential mobility analyzer response functions"</i> (Version 2.0), [Software], Zenodo, https://doi.org/10.5281/zenodo.2652893.

[travis_badge]: https://travis-ci.org/mdpetters/DifferentialMobilityAnalyzers.jl.svg?branch=master
[travis_url]: https://travis-ci.org/mdpetters/DifferentialMobilityAnalyzers.jl

[appveyor_badge]: https://ci.appveyor.com/api/projects/status/github/mdpetters/DifferentialMobilityAnalyzers.jl?svg=true&branch=master

[appveyor_url]: https://ci.appveyor.com/project/mdpetters/differentialmobilityanalyzers-jl

[codecov_badge]: http://codecov.io/github/mdpetters/DifferentialMobilityAnalyzers.jl/coverage.svg?branch=master
[codecov_url]: http://codecov.io/github/mdpetters/DifferentialMobilityAnalyzers.jl?branch=master

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://mdpetters.github.io/DifferentialMobilityAnalyzers.jl/latest/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://mdpetters.github.io/DifferentialMobilityAnalyzers.jl/stable/
