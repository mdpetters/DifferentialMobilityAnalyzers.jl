# DifferentialMobilityAnalyzers.jl NEWS

#### Notes on release changes and ongoing development

- current master is compatible with Julia 1.5
- v2.1 supports Julia 1.5
- v2.0 supports Julia 1.1
- Version 1.0.0 is the last one to support Julia 0.6.4

---
## (current master)

## 2.1
- Fix deprecation warnings from DataFrames API
- Use generic types: AbstractFloat, Vector{<:AbstractFloat} and AbstractMatrix. This fixes an error for  𝐀 * 𝕟 if 𝐀 is of Adjoint type, the default output of Eq. 8 in the paper.
- Ensure symmetry to size distribution arithmetic, i.e. 𝐀 * 𝕟 == 𝕟 * 𝐀; a * 𝕟 == 𝕟 * a; T * 𝕟 == 𝕟 * T; a ⋅ 𝕟 == 𝕟 ⋅ a. 
- Update packages and tests to be compatible with Julia 1.5

## 2.0
- Fix deprecations to be compatible with Julia 1 series.
- Switch all Jupyter documentation graphs from javascript to svg. 
- Remove broadcasting . operator for all operations with size distributions (𝕟 * 𝕟 instead of 𝕟 .* 𝕟 and EF ⋅ 𝕟 instead of EF .⋅ 𝕟). 
- Designate DMA geometry :cylindrical or :radial in DMAconfig   

## 1.0.0
- Version released with the manuscript
