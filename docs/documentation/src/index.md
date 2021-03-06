# DifferentialMobilityAnalyzers.jl

*A Julia package for working with data from differential mobility analyzers.*

DifferentialMobilityAnalyzers.jl bundles a set of abstractions to write concise forward and inverse models of experimental setups that involve  differential mobility analyzers. Abstractions include specialized data types, operators, functions, and conventions. 
Conventions include rules for typesetting font and sub- and superscripting variables. 

If you use this package in your research, please cite the relevant sources listed below. 

## Package Features

- Primitives to simplify working with size distribution data
- Implementation of physical equations describing DMAs
- Computational solution of the discretized Fredholm integral equation
- Generation of custom convolution matrices 
- Fast size distribution inversion using Tikhonov regularization
- Modeling of particle transmission through single and chained DMA setups

## Citations

This work was supported by the United States Department of Energy, Office of Science, Biological and Environment Research, Grant number DE-SC0018265. The tutorial was supported by the American Association of Aerosol Research. 

Petters, M.D. (2018) _A language to simplify computation of differential mobility analyzer response functions_ Aerosol Science & Technology, 52 (12), 1437-1451, [https://doi.org/10.1080/02786826.2018.1530724](https://doi.org/10.1080/02786826.2018.1530724).

Petters, M.D. (2019, April 27) _Virtual Machine containing Software for "A language to simplify computation of differential mobility analyzer response functions"_ (Version 2.0), [Software], Zenodo, [https://doi.org/10.5281/zenodo.2652893](https://doi.org/10.5281/zenodo.2652893).