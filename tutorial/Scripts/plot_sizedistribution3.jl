using DataFrames

function getplots3(𝕟, 𝕩)
    df1 = DataFrame(Dp = 𝕟.Dp, S = 𝕟.S, Dist = ["𝕟" for i = 1:length(𝕟.Dp)])
    df2 = DataFrame(Dp = 𝕩.Dp, S = 𝕩.S, Dist = ["𝕩" for i = 1:length(𝕩.Dp)])
    df = [df1; df2]

    xlabels = log10.([40, 100, 400])
    colors = ["black", "steelblue3"]

    set_default_plot_size(16cm, 9cm)
    return plot(
        df,
        x = :Dp,
        y = :S,
        color = :Dist,
        Geom.step,
        Guide.xlabel("Particle diameter (nm)"),
        Guide.ylabel("dN/dlnD (cm-3)"),
        Guide.xticks(ticks = log10.([40, 50, 60, 70, 80, 90, 100, 200, 300, 400])),
        Guide.colorkey(; title = "Distribution"),
        Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
        Scale.color_discrete_manual(colors...),
        Coord.cartesian(xmin = log10(40), xmax = log10(400)),
    )
end

getplots3(𝕟, 𝕩)
