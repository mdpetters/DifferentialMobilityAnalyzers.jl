
using DataFrames,
    Gadfly, CSV, DifferentialMobilityAnalyzers, Printf, LsqFit, SpecialFunctions

t, p, lpm = 293.15, 940e2, 1.666e-5
r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369
Λ = DMAconfig(t, p, 1lpm, 4lpm, r₁, r₂, l, 0.0, :+, 6, :cylindrical)
δ = setupDMA(Λ, vtoz(Λ, 10000), vtoz(Λ, 10), 120)

# Load a simple comma delimited text file - file contains :Dp, :Rcn, :Rccn
# Note that this is the same DMA/size distribution defined by this DMA
df = CSV.read("example_data.csv")
𝕣ᶜⁿ = (df, :Dp, :Rcn, δ) |> interpolateDataFrameOntoδ        # CN response distribution
𝕣ᶜᶜⁿ = (df, :Dp, :Rccn, δ) |> interpolateDataFrameOntoδ;     # CCN response distribution

function threshold!(𝕟::SizeDistribution, c::Float64, n1::Float64, n2::Float64)
    N = 𝕟.N
    S = 𝕟.S
    S[(N.<=c).&(𝕟.Dp.>150)] .= n2
    N[(N.<=c).&(𝕟.Dp.>150)] .= n2
    S[(N.<=c).&(𝕟.Dp.<150)] .= n1
    N[(N.<=c).&(𝕟.Dp.<150)] .= n1
    𝕟.N = N
end

threshold!(𝕣ᶜⁿ, 0.1, 0.1, 0.1)
threshold!(𝕣ᶜᶜⁿ, 0.1, 0.0, 0.1)

𝕒𝕗 = 𝕣ᶜᶜⁿ / 𝕣ᶜⁿ

Taf(𝕟, μ, σ) = @. 0.5 * (1.0 + erf((𝕟.Dp - μ) ./ (sqrt(2.0σ))));

𝐈, 𝐒, 𝐀, λ = δ.𝐈, δ.𝐒, δ.𝐀, 0.5
𝕟ᶜⁿ = (𝐀'𝐀 + λ^2𝐈)^(-1) * (𝐀'𝕣ᶜⁿ + λ^2 * 𝐒^(-1) * 𝕣ᶜⁿ)
model(x, p) = (𝐀 * (𝕟ᶜⁿ.N .* Taf(𝕟ᶜⁿ, p[1], p[2]))) ./ (𝐀 * 𝕟ᶜⁿ.N)
fit = curve_fit(model, 𝕒𝕗.Dp, 𝕒𝕗.N, [65.0, 3.0])
Ax = fit.param
afmodel = model(δ.Dp, Ax)

df1 = DataFrame(Dp = 𝕣ᶜⁿ.Dp, S = 𝕣ᶜⁿ.S, Dist = ["𝕣ᶜⁿ" for i = 1:length(𝕣ᶜⁿ.Dp)])
df2 = DataFrame(Dp = 𝕣ᶜᶜⁿ.Dp, S = 𝕣ᶜᶜⁿ.S, Dist = ["𝕣ᶜᶜⁿ" for i = 1:length(𝕣ᶜᶜⁿ.Dp)])
df = [df1; df2]

dfr1 = DataFrame(Dp = 𝕒𝕗.Dp, S = 𝕒𝕗.N, Dist = ["𝕒𝕗 (data)" for i = 1:length(𝕒𝕗.Dp)])
dfr2 = DataFrame(Dp = 𝕒𝕗.Dp, S = afmodel, Dist = ["𝕒𝕗 (model)" for i = 1:length(𝕒𝕗.Dp)])
dfr = [dfr1; dfr2]

xlabels = log10.([10, 100, 500])
colors = ["darkred", "steelblue3", "black"]
p1 = plot(
    dfr,
    x = :Dp,
    y = :S,
    color = :Dist,
    Geom.step,
    Guide.xlabel("Apparent +1 Mobility Diameter (nm)"),
    Guide.ylabel("Fraction (-)"),
    Guide.title("Activated Fraction"),
    Guide.colorkey(; title = ""),
    Guide.xticks(
        ticks = log10.([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500]),
    ),
    Guide.yticks(ticks = collect(0:0.2:1.2)),
    Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
    Scale.color_discrete_manual(["black", "darkgoldenrod"]...),
    Coord.cartesian(xmin = log10(10), xmax = log10(500), ymin = 0, ymax = 1.2),
)

p2 = plot(
    df,
    x = :Dp,
    y = :S,
    color = :Dist,
    Geom.step,
    Guide.xlabel("Apparent +1 Mobility Diameter (nm)"),
    Guide.ylabel("dN/dlnD (cm-3)"),
    Guide.title("Raw response function"),
    Guide.xticks(
        ticks = log10.([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500]),
    ),
    Guide.colorkey(; title = ""),
    Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
    Scale.color_discrete_manual(colors...),
    Coord.cartesian(xmin = log10(10), xmax = log10(500), ymin = 0),
)
set_default_plot_size(20cm, 8cm)
p = (hstack(p2, p1))
