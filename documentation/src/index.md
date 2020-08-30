# DifferentialMobilityAnalyzers.jl

*A Julia package for working with data from differential mobility analyzers.*

DifferentialMobilityAnalyzers.jl bundles a set of abstractions to write concise forward and inverse models of experimental setups that involve  differential mobility analyzers. Abstractions include specialized data types, operators, functions, and conventions. Conventions include rules for typesetting fonts and sub- and superscripting variables. 

## Package Features

- Primitives to simplify working with size distribution data
- Implementation of physical equations describing DMAs
- Computational solution of the discretized Fredholm integral equation
- Generation of custom convolution matrices 
- Fast size distribution inversion using Tikhonov regularization
- Modeling of particle transmission through single and chained DMA setups