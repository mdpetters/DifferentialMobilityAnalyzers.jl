# DifferentialMobilityAnalyzers.jl NEWS

#### Notes on release changes and ongoing development

- (current master) supports Julia 1.1
- Version 1.0.0 is the last one to support Julia 0.6.4

---
## (current master)
- Fix deprecations for Julia 1 series
- Switch all Jupyter documentation graphs from javascript to svg. 
- Remove broadcasting . operator for all operations with size distributions (𝕟 * 𝕟 instead of 𝕟 .* 𝕟 and EF ⋅ 𝕟 instead of EF .⋅ 𝕟). 
- Designate DMA geometry :cylindrical or :radial in DMAconfig   

## 1.0.0
- Version released with the manuscript