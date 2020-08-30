using Gadfly, Printf

function getTc_plot()
    t, p = 295.15, 1e5
    qsa, β = 1.66e-5, 1 / 5
    r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369
    leff = 13.0
    m = 6
    DMAtype = :cylindrical
    polarity = :-

    Λp = DMAconfig(t, p, qsa, qsa / β, r₁, r₂, l, leff, :-, m, DMAtype)
    bins, z₁, z₂ = 2, vtoz(Λp, 10000), vtoz(Λp, 10)
    δp = setupDMA(Λp, z₁, z₂, bins)

    Λm = DMAconfig(t, p, qsa, qsa / β, r₁, r₂, l, leff, :+, m, DMAtype)
    bins, z₁, z₂ = 2, vtoz(Λm, 10000), vtoz(Λm, 10)
    δm = setupDMA(Λm, z₁, z₂, bins)


    Dp = exp10.(range(log10(1.0), stop = log10(1000.0), length = 100))
    δp.Tc(1, Dp)

    df1 = DataFrame(Dp = Dp, Tc = δp.Tc(1, Dp), k = ["+1" for i = 1:length(Dp)])
    df2 = DataFrame(Dp = Dp, Tc = δm.Tc(1, Dp), k = ["-1" for i = 1:length(Dp)])
    df3 = DataFrame(Dp = Dp, Tc = δp.Tc(2, Dp), k = ["+2" for i = 1:length(Dp)])
    df4 = DataFrame(Dp = Dp, Tc = δm.Tc(2, Dp), k = ["-2" for i = 1:length(Dp)])
    df5 = DataFrame(Dp = Dp, Tc = δm.Tc(3, Dp), k = ["±3" for i = 1:length(Dp)])
    df6 = DataFrame(Dp = Dp, Tc = δm.Tc(4, Dp), k = ["±4" for i = 1:length(Dp)])
    df = [df2; df1; df4; df3; df5; df6]

    gengrid(r) = [vcat(map(x -> x:x:9x, r)...); r[end] * 10]
    xlabels = log10.([10, 20, 50, 100, 200, 500, 1000])
    ylabels = log10.([0.02, 0.05, 0.1, 0.2, 0.5])
    colors = ["black", "grey", "darkred", "salmon1", "steelblue3", "darkgoldenrod"]
    set_default_plot_size(16cm, 9cm)

    return plot(
        df,
        x = :Dp,
        y = :Tc,
        color = :k,
        Geom.line,
        Guide.xlabel("Particle diameter (nm)"),
        Guide.ylabel("Tc (-)"),
        Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
        Scale.y_log10(labels = y -> y in ylabels ? @sprintf("%.2f", exp10(y)) : ""),
        Guide.xticks(ticks = log10.(gengrid([10, 100]))),
        Guide.yticks(ticks = log10.([0.02:0.01:0.09; 0.1:0.1:0.5])),
        Scale.color_discrete_manual(colors...),
        Coord.cartesian(
            xmin = log10(10),
            xmax = log10(1000),
            ymin = log10(0.02),
            ymax = log10(0.5),
        ),
    )
end

getTc_plot()
