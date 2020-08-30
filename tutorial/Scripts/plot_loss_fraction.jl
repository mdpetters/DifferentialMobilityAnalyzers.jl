
function loss_fraction()
    t, p = 295.15, 1e5
    r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369
    leff = 13.0
    m = 6
    DMAtype = :cylindrical
    polarity = :-
    lpm = 16.66666 * 1e-6

    Dp = exp10.(range(log10(1.0), stop = log10(1000.0), length = 100))

    Λp = DMAconfig(t, p, 1lpm, 10lpm, r₁, r₂, l, leff, :-, m, DMAtype)
    bins, z₁, z₂ = 2, vtoz(Λp, 10000), vtoz(Λp, 10)
    δp = setupDMA(Λp, z₁, z₂, bins)
    df1 =
        DataFrame(Dp = Dp, Tl = δp.Tl(Λp, Dp), leff = ["1:10 L min-1" for i = 1:length(Dp)])

    Λp = DMAconfig(t, p, 0.3lpm, 3lpm, r₁, r₂, l, leff, :-, m, DMAtype)
    bins, z₁, z₂ = 2, vtoz(Λp, 10000), vtoz(Λp, 10)
    δp = setupDMA(Λp, z₁, z₂, bins)
    df2 = DataFrame(
        Dp = Dp,
        Tl = δp.Tl(Λp, Dp),
        leff = ["0.3:3 L min-1" for i = 1:length(Dp)],
    )

    df = [df1; df2]
    gengrid(r) = [vcat(map(x -> x:x:9x, r)...); r[end] * 10]
    xlabels = log10.([1, 10, 100, 1000])
    ylabels = log10.([0.01, 0.1, 1])
    colors = ["black", "darkred"]
    set_default_plot_size(16cm, 9cm)


    plot(
        df,
        x = :Dp,
        y = :Tl,
        color = :leff,
        Geom.line,
        Guide.xlabel("Particle diameter (nm)"),
        Guide.ylabel("Tl (-)"),
        Guide.xticks(ticks = log10.(gengrid([1, 10, 100]))),
        Guide.yticks(ticks = log10.(gengrid([0.01, 0.1]))),
        Guide.colorkey(; title = "Flow rates"),
        Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
        Scale.y_log10(labels = y -> y in ylabels ? @sprintf("%.2f", exp10(y)) : ""),
        Scale.color_discrete_manual(colors...),
        Coord.cartesian(
            xmin = log10(1),
            xmax = log10(1000),
            ymin = log10(0.01),
            ymax = log10(1),
        ),
    )
end

loss_fraction()
