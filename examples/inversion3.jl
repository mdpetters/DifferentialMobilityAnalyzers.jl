using DifferentialMobilityAnalyzers
using Gadfly
using LinearAlgebra
using Printf
using Distributions
using DataFrames
using MLStyle
using CSV
using Lazy
using Colors
using Underscores
using Random
using RegularizationTools
	
function plot_gf_dual(dfl, dfr)
    colors = ["black", "darkred", "steelblue3", "darkgoldenrod"]
    xlabels = collect(1:0.5:3)
    p1 = plot(dfl, x = :gf, y = :N, color = :Color, Geom.step,
        Theme(plot_padding=[0mm, 10mm, 2mm, 2mm]), 
        Guide.xlabel("Apparent Growth Factor (-)"),
        Guide.ylabel("Concentration (cm⁻³)", orientation = :vertical),
        Guide.xticks(ticks = collect(0.8:0.1:2.5)),
        Guide.colorkey(title = ""),
        Scale.color_discrete_manual(colors...),
        Scale.x_continuous(labels = x -> x in xlabels ? @sprintf("%.1f", x) : ""),
        Coord.cartesian(xmin = 0.8, xmax = 2.5))

    p2 = plot(dfr, x = :gf, y = :Frequency, color = :Color, Geom.step,
        Theme(plot_padding=[-8mm, 2mm, 2mm, 2mm]), 
        Guide.xlabel("Growth Factor (-)"),
        Guide.ylabel("Probability Density (-)", orientation = :vertical),
        Guide.xticks(ticks = collect(0.8:0.1:2.5)),
        Guide.colorkey(title = ""),
        Scale.color_discrete_manual(colors...),
        Scale.x_continuous(labels = x -> x in xlabels ? @sprintf("%.1f", x) : ""),
        Coord.cartesian(xmin = 0.8, xmax = 2.5))
	
    hstack(p1,p2)
end

include("commonTDMAfunctions.jl")

Qcpc = 1.0 # Flow in LPM
Dd = 100e-9
Nt = 3000.0
k = 60
seed = 2000
gf0 = 1.6           

Λ₁, Λ₂, δ₁, δ₂ = initializeDMAs(Dd, k)
Ax = [[1300.0, 60.0, 1.4], [2000.0, 200.0, 1.6]]
𝕟ᶜⁿ = DMALognormalDistribution(Ax, δ₁)
gf, ge, 𝐀 = TDMAmatrix(𝕟ᶜⁿ, Dd, Λ₁, Λ₂, δ₂, k)
model = TDMA1Dpdf(𝕟ᶜⁿ, Λ₁, Λ₂, (Dd, 0.8, 5.0, k));

dg = ge[1:end-1] .- ge[2:end]
f = @> (0.7*pdf.(Normal(1.3,0.1), gf) + pdf.(Normal(1.7,0.25), gf)) Normalize
N0 = 𝐀*f
N1 = poisson_noise(Qcpc, N0; seed = seed, t = 2.0)

x₀ = N1./sum(N1)
lb, ub = zeros(k), ones(k)
xλ1 = invert(𝐀, N1, Lₖx₀B(2, x₀, lb, ub))
e1 = @> sqrt.(sum((xλ1 .- f).^2.0)./k) round(digits = 3)
xλ2 = invert(𝐀, N1, LₖDₓB(0, 0.001, lb, ub))
e2 = @> sqrt.(sum((xλ2 .- f).^2.0)./k) round(digits = 3)
 
dfl2 = DataFrame(gf = gf, Frequency = xλ2 ./ dg, Color = "L<sub>0</sub>D<sub>1e-3</sub>B<sub>[0,1]</sub>, $(e2)")
dfl3 = DataFrame(gf = gf, Frequency = xλ1 ./ dg, Color = "L<sub>2</sub>x<sub>0</sub>B<sub>[0,1]</sub>, $(e1)")

dfr1 = DataFrame(gf = gf, N = N0, Color = "No noise")
dfr2 = DataFrame(gf = gf, N = N1, Color = "Q = 1lpm")
dfl1 = DataFrame(gf = gf, Frequency = f ./ dg, Color = "Truth")

set_default_plot_size(18cm, 7cm)
p = plot_gf_dual([dfr1;dfr2], [dfl1; dfl2; dfl3])
#Gadfly.draw(SVG("f03.svg"), p)
