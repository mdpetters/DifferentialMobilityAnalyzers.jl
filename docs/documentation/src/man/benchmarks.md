# Performance Benchmarks

Currently best performance is achieved running on a single BLAS threads (set as default) and the using the [Intel MKL linear algebra library](https://github.com/JuliaComputing/MKL.jl). The MKL library is not shipped in the default Julia distribution must be manually installed.

 [Benchmarks](@ref) can be computed using supplied functions. Running the benchmark test will take some time since the routine takes an average of multiple runs and cycles through several configurations.

The three slowest routines included in the benchmark are [rinv](@ref), [setupDMA](@ref), and [setupSMPS](@ref). Setting up the DMAs is a one-time computational cost. If speed is absolutely critical, the fields Λ and δ can be precomputed and loaded from file. The computational cost of [rinv](@ref) is recurring.

Example benchmark run:

```@example
using DifferentialMobilityAnalyzers
runbenchmarks()
```