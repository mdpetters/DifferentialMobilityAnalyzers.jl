using Distributions
using DifferentialMobilityAnalyzers
using Gadfly

t, p = 295.15, 1e5
qsa, qsh = 1.66e-5, 8.33e-5
r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369
Λ₁ = DMAconfig(t, p, qsa, qsh, r₁, r₂, l, 0.0, :-, 6, :cylindrical)
Λ₂ = DMAconfig(t, p, qsa, qsh, r₁, r₂, l, 0.0, :-, 6, :cylindrical)
bins, z₁, z₂ = 120, dtoz(Λ₁, 500e-9), dtoz(Λ₁, 30e-9) # bins, upper, lower mobility limit
δ₁ = setupDMA(Λ₁, z₁, z₂, bins)                

Ax = [[1300.0, 60.0, 1.4], [5000.0, 220.0, 1.6]] 
𝕟 = DMALognormalDistribution(Ax, δ₁)

# scan 100 nm Dd from 0.8Dd to 3.0Dd with 100 bins
dma2range = (100e-9, 0.8, 3.0, 120)

# Get the model function
model = TDMA1Dpdf(𝕟, Λ₁, Λ₂, dma2range)

P = [0.5,0.15, 0.10, 0.25]   # Probability of growth factor (4 populations)
gf = [1.0, 1.2, 1.6, 2.1]    # Values of growth factor
𝕘 = model(𝕟, P, dma2range[1], gf)

plot(x = 𝕘.Dp/(dma2range[1]*1e9), y = 𝕘.N, Geom.line,
    Guide.xticks(ticks = 0.8:0.2:3),
    Coord.cartesian(xmin = 0.8, xmax = 3.0))