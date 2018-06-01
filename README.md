## DifferentialMobilityAnalyzers.jl
[![travis badge][travis_badge]][travis_url]
[![appveyor badge][appveyor_badge]][appveyor_url]
[![codecov badge][codecov_badge]][codecov_url]
[![][docs-stable-img]](docs/DOCUMENTATION.md)

This package implements the Julia DMA language. The language is a tool to simplify interpretation of data from aerosol differential mobility analyzers.

## Installation & Quickstart
The package  <b> DifferentialMobiliyAnalyzers </b> can be installed from the Julia REPL prompt with
```julia
julia> Pkg.clone("https://github.com/mdpetters/DifferentialMobilityAnalyzers.jl.git")
```
This installs the package and any missing dependencies. (Patience required for fresh install).
```julia
julia> using IJulia
julia> notebook(detached = true)
```
This starts the notebook server. Then load any of the notebooks in the docs/ folder.

## Documentation
Detailed documentation: [![][docs-stable-img]](docs/DOCUMENTATION.md)

Quickstart: The Julia DMA language is documented in a manuscript and 12 Supplementary Jupyter Notebooks. The links open the notebooks in viewer mode via NBViewer

[Julia DMA Language (preprint)]()<br>
[Notebook S1. Differential Mobility Analyzer](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S01.%20Differential%20Mobility%20Analyzer.ipynb) <br>
[Notebook S2. Fredholm Integral Equation](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S02.%20Fredholm%20Integral%20Equation.ipynb) <br>
[Notebook S3. Size Distribution Arithmetic](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S03.%20Size%20Distribution%20Arithmetic.ipynb) <br>
[Notebook S4. Single Mobility Classification](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S04.%20Single%20Mobility%20Classification.ipynb) <br>
[Notebook S5. Size Distribution Inversion Using Regularization](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S05.%20Size%20Distribution%20Inversion%20Using%20Regularization.ipynb) <br>
[Notebook S6. Size Distribution Inversion of Ambient Data](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S06.%20Size%20Distribution%20Inversion%20of%20Ambient%20Data.ipynb) <br>
[Notebook S7. Size resolved CCN measurements](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S07.%20Size%20resolved%20CCN%20measurements.ipynb) <br>
[Notebook S8. Hygroscopicity Tandem DMA](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S08.%20Hygroscopicity%20Tandem%20DMA.ipynb) <br>
[Notebook S9. Volatility Tandem DMA](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S09.%20Volatility%20Tandem%20DMA.ipynb) <br>
[Notebook S10. Dimer Coagulation and Isolation](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S10.%20Dimer%20Coagulation%20and%20Isolation.ipynb) <br>
[Notebook S11. PartMC Simulations](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S11.%20PartMC%20Simulations.ipynb)<br>
[Notebook S12. FORTRAN API](https://nbviewer.jupyter.org/github/mdpetters/DifferentialMobilityAnalyzers.jl/blob/master/docs/Notebook%20S12.%20FORTRAN%20API.ipynb) <br>

## Contribute
Contributions including notebooks for classroom instruction, homework assignments, interesting DMA configurations, new inversion schemes, and improved or new functionalities of the language are welcome.

## Citation
This work was supported by the United States Department of Energy, Office of Science, Biological and Environment Research, Grant number DE-SC0018265. If you use <b> DifferentialMobilityAnalyzers.jl </b> in your research, please cite

Petters, M.D. (2018) <i> A language to simplify computation of differential mobility analyzer response functions </i> submitted to Aerosol Science & Technology.

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg

[travis_badge]: https://travis-ci.org/mdpetters/DifferentialMobilityAnalyzers.jl.svg?branch=master
[travis_url]: https://travis-ci.org/mdpetters/DifferentialMobilityAnalyzers.jl

[appveyor_badge]: https://ci.appveyor.com/api/projects/status/github/mdpetters/DifferentialMobilityAnalyzers.jl?svg=true&branch=master

[appveyor_url]: https://ci.appveyor.com/project/mdpetters/differentialmobilityanalyzers-jl

[codecov_badge]: http://codecov.io/github/mdpetters/DifferentialMobilityAnalyzers.jl/coverage.svg?branch=master
[codecov_url]: http://codecov.io/github/mdpetters/DifferentialMobilityAnalyzers.jl?branch=master
