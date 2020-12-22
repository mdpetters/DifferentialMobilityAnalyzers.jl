# Transmission Through the Humified tandem DMA

using DifferentialMobilityAnalyzers
using Gadfly
using NumericIO
using Colors
using LinearAlgebra
using Printf
using DataFrames
using Underscores
import Lazy.@>, Lazy.@>>

t, p = 295.15, 1e5                                # Temperature [K], Pressure [Pa]
qsa, β = 1.66e-5, 1 / 5                             # Qsample [m3 s-1], Sample-to-sheath ratio,
r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369               # DMA geometry [m]
Λ₁ = DMAconfig(t, p, qsa, qsa / β, r₁, r₂, l, 0.0, :-, 3, :cylindrical)  # Specify DMA1
Λ₂ = DMAconfig(t, p, qsa, qsa / β, r₁, r₂, l, 0.0, :-, 3, :cylindrical)  # Specify DMA2
bins, z₁, z₂ = 512, dtoz(Λ₁, 500e-9), dtoz(Λ₁, 30e-9) # bins, upper, lower mobility limit
δ₁ = setupDMA(Λ₁, z₁, z₂, bins)                  # Compute matrices
δ₂ = setupDMA(Λ₂, z₁, z₂, bins)                  # Compute matrices

# Upstream Size Distribution
Ax = [[1300.0, 60.0, 1.4], [5000.0, 220.0, 1.6]]
𝕟ᶜⁿ = DMALognormalDistribution(Ax, δ₁)

# Tandem DMA equations
T₁(zˢ, k) = δ₁.Ω(Λ₁, δ₁.Z, zˢ / k) .* δ₁.Tc(k, δ₁.Dp) .* δ₁.Tl(Λ₁, δ₁.Dp)
cr(zˢ, k) = ztod(Λ₁, 1, zˢ) / ztod(Λ₁, k, zˢ)
DMA₁(𝕟, zˢ, gf) = @_ map(cr(zˢ, _) ⋅ (gfₖ(Λ₁, zˢ, gf, _) ⋅ (T₁(zˢ, _) * 𝕟)), 1:3)
itp(𝕟) = interpolateSizeDistributionOntoδ((𝕟, δ₂))
DMA₂(𝕟) = δ₂.𝐎 * 𝕟

Dd = 100e-9             # Dry diameter
zˢ = dtoz(Λ₁, Dd);      # Mobility of 100 nm particle
gf = 1.6                # Growth factor
𝕄 = @_ map(itp(_) |> DMA₂, DMA₁(𝕟ᶜⁿ, zˢ, gf)) # 𝕄[k] distributions
𝕞ᵗ = sum(𝕄)                                  # total response

mdf(k) = DataFrame(
    Dp = 𝕄[k].Dp./(Dd*1e9), 
    S = 𝕄[k].S, 
    Dist = ["𝕄[$k]" for i = 1:length(𝕄[k].Dp)]
)

df1 = mapreduce(mdf, vcat, 1:3)
df2 = DataFrame(Dp = 𝕞ᵗ.Dp./(Dd*1e9), S = 𝕞ᵗ.S, Dist = ["𝕞ᵗ" for i = 1:length(𝕞ᵗ.Dp)])
df = [df2; df1]

colors = ["black", "darkred", "steelblue3", "darkgoldenrod"]

p2 = plot(
    df,
    x = :Dp,
    y = :S,
    color = :Dist,
    Geom.line,
    Guide.xlabel("Apparent Growth Factor", orientation = :horizontal),
    Guide.ylabel("dN/dlnD (cm⁻³)"),
    Guide.xticks(ticks = [1.2,1.4,1.6,1.8,2.0]),
    Guide.colorkey(; title = ""),
    Scale.color_discrete_manual(colors...),
    Coord.cartesian(xmin = 1.2, xmax = 2),
    Theme(plot_padding = [5mm, 10mm, 0mm, 0mm]),
)
