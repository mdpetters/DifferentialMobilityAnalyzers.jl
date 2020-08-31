# Transmission Through the DMA at a Constant Voltage

using DifferentialMobilityAnalyzers
using Gadfly
using NumericIO
using Colors
using LinearAlgebra
using Printf
using DataFrames

# Create a DMA config
qsa,qsh = 1.66e-5, 8.33e-5                       # Qsample [m3 s-1], Qsheath [m3 s-1]
t,p = 295.15, 1e5                                # Temperature [K], Pressure [Pa]
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369               # DMA geometry [m]
Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,0.0,:-,6,:cylindrical)

# Create a DMA grid
z₁,z₂ = vtoz(Λ,10000), vtoz(Λ,10)    # bins, upper, lower mobility limit
δ = setupDMA(Λ, z₁, z₂, 512);

# Compute the transmission through the DMA
T(zˢ,k,Λ,δ) = δ.Ω(Λ,δ.Z,zˢ/k).*δ.Tc(k,δ.Dp).*δ.Tl(Λ,δ.Dp)
zˢ = dtoz(Λ, 100*1e-9)
𝕟ᶜⁿ = DMALognormalDistribution([[900., 40., 1.5], [500., 180., 1.4]], δ)
ℕ = map(k -> T(zˢ,k,Λ,δ)*𝕟ᶜⁿ,1:3)
𝕄 = map(k -> (ztod(Λ,1,zˢ)/ztod(Λ,k,zˢ))⋅(T(zˢ,k,Λ,δ)*𝕟ᶜⁿ),1:3)
𝕟ₜ, 𝕞ₜ = sum(ℕ), sum(𝕄)

# Plot the results
set_default_plot_size(25cm, 7cm) 

xlabels = log10.([10, 50, 100, 500]) 
p1 = plot(
	x = 𝕟ᶜⁿ.Dp, 
    y = 𝕟ᶜⁿ.S, 
    Geom.step, 
    color = ["𝕟ᶜⁿ" for i in 𝕟ᶜⁿ.Dp], 
    Guide.xlabel("Particle diameter (nm)"), #
    Guide.ylabel("dN/dlnD (cm-3)"), 
    Guide.xticks( # hide
        ticks = log10.([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600]), 
    ), 
    Guide.colorkey(; title = ""),
    Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
    Scale.color_discrete_manual("black"),
    Coord.cartesian(xmin = log10(10), xmax = log10(600)),
    Theme(plot_padding = [0mm, 0mm, 0mm, 0mm]),
) 

df1 = let
    xx = map(1:3) do i
        df = DataFrame(x = 𝕄[i].Dp, y = 𝕄[i].S, c = ["𝕄[$i]" for j in 𝕄[i].Dp])
    end
    vcat(xx...)
end

df2 = DataFrame(x = 𝕞ₜ.Dp, y = 𝕞ₜ.S, c = ["𝕞ₜ" for j in 𝕞ₜ.Dp])

df = [df2; df1]
colors = ["black", "darkred", "steelblue3", "darkgoldenrod"]# hide

p2 = plot(
    df,
    x = :x,
    y = :y,
    color = :c,
    Geom.line,
    Guide.xlabel("Apparent +1 Mobility Diameter (nm)", orientation = :horizontal),
    Guide.ylabel(""),
    Guide.xticks(ticks = [80, 100, 120, 140]),
    Guide.colorkey(; title = ""),
    Scale.color_discrete_manual(colors...),
    Coord.cartesian(xmin = 70, xmax = 150),
    Theme(plot_padding = [5mm, 10mm, 0mm, 0mm]),
)

df1 = let
    xx = map(1:3) do i
        df = DataFrame(x = ℕ[i].Dp, y = ℕ[i].S, c = ["ℕ[$i]" for j in ℕ[i].Dp])
    end
    vcat(xx...)
end

df2 = DataFrame(x = 𝕟ₜ.Dp, y = 𝕟ₜ.S, c = ["𝕟ₜ" for j in 𝕟ₜ.Dp])

df = [df2; df1]

p3 = plot(
    df,
    x = :x,
    y = :y,
    color = :c,
    Geom.line,
    Guide.xlabel("Mobility Diameter (nm)", orientation = :horizontal),
    Guide.ylabel(""),
    Guide.xticks(ticks = [100, 150, 200, 250]),
    Guide.colorkey(; title = ""),
    Scale.color_discrete_manual(colors...),
    Coord.cartesian(xmin = 70, xmax = 250),
    Theme(plot_padding = [0mm, 10mm, 0mm, 0mm]),
)

hstack(p1, p2, p3)
