
# # Size distribution
using NumericIO
set_default_plot_size(25cm, 7cm)

xlabels = log10.([10, 50, 100, 500])
p1 = plot(
    x = 𝕟ᶜⁿ.Dp,
    y = 𝕟ᶜⁿ.S,
    Geom.step,
    color = ["𝕟ᶜⁿ" for i in 𝕟ᶜⁿ.Dp],
    Guide.xlabel("Particle diameter (nm)"),
    Guide.ylabel("dN/dlnD (cm-3)"),
    Guide.xticks(
        ticks = log10.([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600]),
    ),
    Guide.colorkey(; title = ""),
    Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
    Scale.color_discrete_manual("black"),
    Coord.cartesian(xmin = log10(10), xmax = log10(600)),
    Theme(plot_padding = [0mm, 0mm, 0mm, 0mm]),
)

df1 = let
    xx = map(1:3) do i
        df = DataFrame(x = 𝕄[i].Dp, y = 𝕄[i].S, c = ["𝕄[$i]" for j in 𝕄[i].Dp])
    end
    vcat(xx...)
end

df2 = DataFrame(x = 𝕞ₜ.Dp, y = 𝕞ₜ.S, c = ["𝕞ₜ" for j in 𝕞ₜ.Dp])

df = [df2; df1]

colors = ["black", "darkred", "steelblue3", "darkgoldenrod"]

p2 = plot(
    df,
    x = :x,
    y = :y,
    color = :c,
    Geom.line,
    Guide.xlabel("Apparent +1 Mobility Diameter (nm)"),
    Guide.ylabel(""),
    Guide.xticks(ticks = [80, 100, 120, 140]),
    Guide.colorkey(; title = ""),
    Scale.color_discrete_manual(colors...),
    Coord.cartesian(xmin = 70, xmax = 150),
    Theme(plot_padding = [0mm, 10mm, 0mm, 0mm]),
)

df1 = let
    xx = map(1:3) do i
        df = DataFrame(x = ℕ[i].Dp, y = ℕ[i].S, c = ["ℕ[$i]" for j in ℕ[i].Dp])
    end
    vcat(xx...)
end

df2 = DataFrame(x = 𝕟ₜ.Dp, y = 𝕟ₜ.S, c = ["𝕟ₜ" for j in 𝕟ₜ.Dp])

df = [df2; df1]

p3 = plot(
    df,
    x = :x,
    y = :y,
    color = :c,
    Geom.line,
    Guide.xlabel("Mobility Diameter (nm)"),
    Guide.ylabel(""),
    Guide.xticks(ticks = [100, 150, 200, 250]),
    Guide.colorkey(; title = ""),
    Scale.color_discrete_manual(colors...),
    Coord.cartesian(xmin = 70, xmax = 250),
    Theme(plot_padding = [0mm, 10mm, 0mm, 0mm]),
)

hstack(p1, p2, p3)
