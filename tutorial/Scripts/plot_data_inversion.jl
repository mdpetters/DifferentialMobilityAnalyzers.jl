using DataFrames

function getplotsinv(𝕣, 𝕟ⁱⁿᵛ¹, 𝕟ⁱⁿᵛ², 𝕟ⁱⁿᵛ³)

    df1 = DataFrame(Dp = 𝕟ⁱⁿᵛ¹.Dp, S = 𝕟ⁱⁿᵛ¹.S, Dist = ["𝕟ⁱⁿᵛ¹" for i = 1:length(𝕟ⁱⁿᵛ¹.Dp)])
    df2 = DataFrame(Dp = 𝕟ⁱⁿᵛ².Dp, S = 𝕟ⁱⁿᵛ².S, Dist = ["𝕟ⁱⁿᵛ²" for i = 1:length(𝕟ⁱⁿᵛ².Dp)])
    df3 = DataFrame(Dp = 𝕟ⁱⁿᵛ³.Dp, S = 𝕟ⁱⁿᵛ³.S, Dist = ["𝕟ⁱⁿᵛ³" for i = 1:length(𝕟ⁱⁿᵛ³.Dp)])
    df = [df1; df2; df3]

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
    p = (hstack(p1, p2))
    set_default_plot_size(20cm, 8cm)
    return p
end

getplotsinv(𝕣, 𝕟ⁱⁿᵛ¹, 𝕟ⁱⁿᵛ², 𝕟ⁱⁿᵛ³)
