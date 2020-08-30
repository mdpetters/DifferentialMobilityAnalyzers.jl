using Gadfly, DifferentialMobilityAnalyzers, DataFrames, Colors

function getΩ1()
    t, p = 295.15, 1e5
    qsa, β = 1.66e-5, 1 / 5
    r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369
    leff = 13.0
    m = 6
    DMAtype = :cylindrical
    polarity = :-

    Λ = DMAconfig(t, p, qsa, qsa / β, r₁, r₂, l, leff, polarity, m, DMAtype)
    bins, z₁, z₂ = 5, vtoz(Λ, 10000), vtoz(Λ, 10)
    δ = setupDMA(Λ, z₁, z₂, bins)

    Z = exp10.(range(log10(1e-9), stop = log10(1e-6), length = 1000))
    z1ˢ = dtoz(Λ, 200e-9)
    Ω1 = δ.Ω(Λ, Z, z1ˢ)
    z2ˢ = dtoz(Λ, 20e-9)
    Ω2 = δ.Ω(Λ, Z, z2ˢ)

    df1 = DataFrame(x = Z ./ z1ˢ, Ω = Ω1, Dp = ["200 nm" for i = 1:length(Z)])
    df2 = DataFrame(x = Z ./ z2ˢ, Ω = Ω2, Dp = ["20 nm" for i = 1:length(Z)])
    df = [df1; df2]
    set_default_plot_size(16cm, 9cm)

    p = plot(
        df,
        x = :x,
        y = :Ω,
        color = :Dp,
        Geom.line,
        Guide.xlabel("z/z <sup>s</sup> (-)"),
        Guide.ylabel("Ω (-)"),
        Guide.xticks(ticks = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]),
        Guide.yticks(ticks = 0:0.2:1),
        Scale.color_discrete_manual(["black", "darkred"]...),
        Coord.cartesian(xmin = 0.7, xmax = 1.3, ymin = 0, ymax = 1),
    )
end

getΩ1()
