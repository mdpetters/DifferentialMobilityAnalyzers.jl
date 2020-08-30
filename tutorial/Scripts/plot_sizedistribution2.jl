using DataFrames

function getplots2(T, 𝕟, 𝕩)
    df1 = DataFrame(Dp = 𝕟.Dp, S = T(𝕟.Dp), Dist = ["T" for i = 1:length(𝕟.Dp)])
    df2 = DataFrame(Dp = 𝕟.Dp, S = 𝕟.S, Dist = ["𝕟" for i = 1:length(𝕟.Dp)])
    df3 = DataFrame(Dp = 𝕩.Dp, S = 𝕩.S, Dist = ["𝕩" for i = 1:length(𝕩.Dp)])
    df = [df2; df3]

    xlabels = log10.([40, 100, 400])
    colors = ["darkred", "steelblue3", "black"]
    p1 = plot(
        df1,
        x = :Dp,
        y = :S,
        color = :Dist,
        Geom.step,
        Guide.xlabel("Particle diameter (nm)"),
        Guide.ylabel("Fraction (-)"),
        Guide.colorkey(; title = ""),
        Guide.xticks(ticks = log10.([40, 50, 60, 70, 80, 90, 100, 200, 300, 400])),
        Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
        Scale.color_discrete_manual("black"),
        Coord.cartesian(xmin = log10(40), xmax = log10(400)),
    )

    p2 = plot(
        df,
        x = :Dp,
        y = :S,
        color = :Dist,
        Geom.step,
        Guide.xlabel("Particle diameter (nm)"),
        Guide.ylabel("dN/dlnD (cm-3)"),
        Guide.xticks(ticks = log10.([40, 50, 60, 70, 80, 90, 100, 200, 300, 400])),
        Guide.colorkey(; title = ""),
        Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
        Scale.color_discrete_manual(colors...),
        Coord.cartesian(xmin = log10(40), xmax = log10(400)),
    )
    p = (hstack(p1, p2))
    set_default_plot_size(20cm, 8cm)
    return p
end

getplots2(T, 𝕟, 𝕩)
