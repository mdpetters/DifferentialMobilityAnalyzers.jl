# Transmission Through the Humidified Tandem DMA 

using DifferentialMobilityAnalyzers
using Gadfly
using NumericIO
using Colors
using LinearAlgebra
using Printf
using DataFrames

t, p = 295.15, 1e5                                # Temperature [K], Pressure [Pa]
qsa, β = 1.66e-5, 1 / 5                             # Qsample [m3 s-1], Sample-to-sheath ratio,
r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369               # DMA geometry [m]
Λ₁ = DMAconfig(t, p, qsa, qsa / β, r₁, r₂, l, 0.0, :-, 3, :cylindrical)  # Specify DMA1
Λ₂ = DMAconfig(t, p, qsa, qsa / β, r₁, r₂, l, 0.0, :-, 3, :cylindrical)  # Specify DMA2
bins, z₁, z₂ = 512, dtoz(Λ₁, 500e-9), dtoz(Λ₁, 30e-9) # bins, upper, lower mobility limit
δ₁ = setupDMA(Λ₁, z₁, z₂, bins)                  # Compute matrices
δ₂ = setupDMA(Λ₂, z₁, z₂, bins)                  # Compute matrices

# Upstream Size Distribution
Ax = [[1300.0, 60.0, 1.4], [2000.0, 200.0, 1.6]]
𝕟ᶜⁿ = DMALognormalDistribution(Ax, δ₁)

# Tandem DMA equations
T(zˢ, k, Λ, δ) = δ.Ω(Λ, δ.Z, zˢ / k) .* δ.Tc(k, δ.Dp) .* δ.Tl(Λ, δ.Dp)
DMA₁(𝕟, zˢ, gf, Λ, δ) =
    sum(map(k -> (ztod(Λ, 1, zˢ) / ztod(Λ, k, zˢ)) ⋅ (gf ⋅ (T(zˢ, k, Λ, δ) * 𝕟)), 1:3))
DMA₂(𝕟, δ) = δ.𝐎 * 𝕟
model(zˢ, gf) = DMA₁(𝕟ᶜⁿ, zˢ, gf, Λ₁, δ₁) |> 𝕟 -> DMA₂(𝕟, δ₂)

# Pass size distribution to DMA₁ and then output distributiom to DMA₂
zˢ = dtoz(Λ₁, 100e-9);   # Mobility of 100 nm particle
gf = 1.55                # Growth factor
𝕞 = model(zˢ, gf)

set_default_plot_size(14cm, 7cm)

xlabels = collect(100:20:240)
p1 = plot(
    x = 𝕞.Dp,
    y = 𝕞.N,
    Geom.step,
    Guide.xlabel("Particle diameter (nm)"),
    Guide.ylabel("Number concentration (cm-3)", orientation = :vertical),
    Guide.xticks(ticks = (collect(100:10:220))),
    Scale.x_continuous(labels = x -> x in xlabels ? @sprintf("%2i", (x)) : ""),
    Scale.color_discrete_manual("black"),
    Coord.cartesian(xmin = 100, xmax = 220),
    Theme(plot_padding = [2mm, 2mm, 2mm, 2mm]),
)
