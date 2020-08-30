using DataFrames

function getplotsccn(𝕣ᶜⁿ, 𝕣ᶜᶜⁿ, 𝕒𝕗)

    p1 = plot(x = 𝕒𝕗.Dp, y = 𝕒𝕗.N, Geom.line)
    p2 = plot(x = 𝕒𝕗.Dp, y = 𝕒𝕗.N, Geom.line)

    df1 = DataFrame(Dp = 𝕣ᶜⁿ.Dp, S = 𝕣ᶜⁿ.S, Dist = ["𝕣ᶜⁿ" for i = 1:length(𝕣ᶜⁿ.Dp)])
    df2 = DataFrame(Dp = 𝕣ᶜᶜⁿ.Dp, S = 𝕣ᶜᶜⁿ.S, Dist = ["𝕣ᶜᶜⁿ" for i = 1:length(𝕣ᶜᶜⁿ.Dp)])
    df = [df1; df2]

    dfr = DataFrame(Dp = 𝕒𝕗.Dp, S = 𝕒𝕗.N, Dist = ["𝕒𝕗" for i = 1:length(𝕣.Dp)])

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
        Scale.color_discrete_manual("black"),
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
    p = (hstack(p2, p1))
    set_default_plot_size(20cm, 8cm)
    return p
end

getplotsccn(𝕣ᶜⁿ, 𝕣ᶜᶜⁿ, 𝕒𝕗)
