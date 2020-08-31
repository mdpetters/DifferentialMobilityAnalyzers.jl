using DataFrames, Gadfly, CSV, DifferentialMobilityAnalyzers, Printf

df = CSV.read("example_data.csv", DataFrame)

t, p, lpm = 293.15, 940e2, 1.666e-5
r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369
Λ = DMAconfig(t, p, 1lpm, 4lpm, r₁, r₂, l, 0.0, :+, 6, :cylindrical)
δ = setupDMA(Λ, vtoz(Λ, 10000), vtoz(Λ, 10), 120);

𝕣 = (df, :Dp, :Rcn, δ) |> interpolateDataFrameOntoδ
𝕟ⁱⁿᵛ = rinv(𝕣.N, δ, λ₁ = 0.1, λ₂ = 1.0)

df = DataFrame(Dp = 𝕟ⁱⁿᵛ.Dp, S = 𝕟ⁱⁿᵛ.S, Dist = ["𝕟ⁱⁿᵛ" for i = 1:length(𝕟ⁱⁿᵛ.Dp)])
dfr = DataFrame(Dp = 𝕣.Dp, S = 𝕣.N, Dist = ["𝕣" for i = 1:length(𝕣.Dp)])

xlabels = log10.([10, 100, 500])
colors = ["darkred", "steelblue3", "black"]
p1 = plot(
    dfr,
    x = :Dp,
    y = :S,
    color = :Dist,
    Geom.step,
    Guide.xlabel("Particle diameter (nm)"),
    Guide.ylabel("N (cm-3)"),
    Guide.title("Raw Response Function"),
    Guide.colorkey(; title = ""),
    Guide.xticks(
        ticks = log10.([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500]),
    ),
    Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
    Scale.color_discrete_manual("black"),
    Coord.cartesian(xmin = log10(10), xmax = log10(500), ymin = 0),
)

p2 = plot(
    df,
    x = :Dp,
    y = :S,
    color = :Dist,
    Geom.step,
    Guide.xlabel("Particle diameter (nm)"),
    Guide.ylabel("dN/dlnD (cm-3)"),
    Guide.title("Inverted Size Distribution"),
    Guide.xticks(
        ticks = log10.([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500]),
    ),
    Guide.colorkey(; title = ""),
    Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
    Scale.color_discrete_manual(colors...),
    Coord.cartesian(xmin = log10(10), xmax = log10(500), ymin = 0),
)
set_default_plot_size(20cm, 8cm)
p = (hstack(p1, p2))
