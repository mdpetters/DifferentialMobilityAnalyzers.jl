## Overview
<b> DifferentialMobilityAnalyzers.jl </b> is an implementation of the Julia DMA language (JDL) and tool for working with data from aerosol differential mobility analyzers.

## Installation & Quickstart
The package  <b> DifferentialMobiliyAnalyzers </b> can be installed from the Julia REPL prompt with
```julia
julia> Pkg.add("DifferentialMobilityAnalyzers")
```
This installs the package and any missing dependencies.
```julia
julia> using IJulia
julia> notebook(detached = true)
```
This starts the notebook server. Then load any of the notebooks in the docs/ folder.

## Documentation
Detailed installation instructions and documentation is [here](docs/DOCUMENTATION.md). The JDL language is described in the manuscript and  12 Supplementary Jupyter Notebooks. The links open the notebooks in viewer mode via NBViewer

[Julia DMA Language (preprint)]()<br>
[Notebook S1. Differential Mobility Analyzer]() <br>
[Notebook S2. Fredholm Integral Equation]() <br>
[Notebook S3. Size Distribution Arithmetic]() <br>
[Notebook S4. Single Mobility Classification]() <br>
[Notebook S5. Size Distribution Inversion Using Regularization]() <br>
[Notebook S6. Size Distribution Inversion of Ambient Data]() <br>
[Notebook S7. Size resolved CCN measurements]() <br>
[Notebook S8. Hygroscopicity Tandem DMA]() <br>
[Notebook S9. Volatility Tandem DMA]() <br>
[Notebook S10. Dimer Coagulation and Isolation]() <br>
[Notebook S11. PartMC Simulations]()<br>
[Notebook S12. FORTRAN API]() <br>

## Contribute
Contributions including notebooks for classroom instruction, homework assignments using JDL, addition of DMA configurations not considered here, new inversion schemes, and improved or new functionalities of the language efficiency or new features) are welcome.

## Citation
This work was supported by the United States Department of Energy, Office of Science, Biological and Environment Research, Grant number DE-SC0018265. If you use <b> DifferentialMobilityAnalyzers.jl </b> in your research, please cite

Petters, M.D. (2018) <i> A language to simplify computation of differential mobility analyzer response functions </i> submitted to Aerosol Science & Technology.
