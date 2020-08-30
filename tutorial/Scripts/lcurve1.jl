function getpsd()
    Random.seed!(703)
    bins = 120
    i =
        (matrices[!, :bins] .== bins) .&
        (matrices[!, :Diff] .== true) .&
        (matrices[!, :Loss] .== true)
    δ = matrices[i, :δ][1]
    𝐀 = δ.𝐀
    𝐀i = inv(δ.𝐀)
    N = 500.0
    𝕟 = DMALognormalDistribution([[0.4 * N, 30, 1.2], [0.5 * N, 110, 1.7]], δ)
    lpm = 16.666666 # 16.66 cm³ s⁻¹ = 1 L min⁻¹
    𝕣 = poisson_noise(1 * lpm, 𝐀 * 𝕟)
    eyeM = Matrix{Float64}(I, length(𝕣.N), length(𝕣.N))
    setupRegularization(δ.𝐀, eyeM, 𝕣.N, inv(δ.𝐒) * 𝕣.N, 1)
    𝕟, 𝕣
end

function lcurvefig()
    n1, r1 = getpsd()
    λopt = lcorner(0.05, 1.0; n = 10, r = 3)
    L1, L2, λs, ii = lcurve(0.01, 2.0; n = 50)
    mL1, mL2 = reginv([0.01, 0.05, 0.1, λopt, 0.5, 1.0, 2.0]; r = :L1L2)
    psd = clean((reginv(λopt, r = :Nλ))[1])

    df1 = DataFrame(Dp = n1.Dp, S = n1.S, Dist = ["𝕟 (true)" for i = 1:length(n1.Dp)])
    rd = round(λopt, digits = 2)
    df2 = DataFrame(
        Dp = n1.Dp,
        S = psd ./ n1.ΔlnD,
        Dist = ["𝕟ⁱⁿᵛ, λopt = $(rd)" for i = 1:length(n1.Dp)],
    )
    dfxx = [df1; df2]

    label = map(x -> @sprintf("λ = %.2f", x), [0.01, 0.05, 0.1, λopt, 0.5, 1.0, 2.0])
    xlabels = log10.([1, 10, 100])
    ylabels = log10.([1, 10])
    colors = ["black", "darkred"]
    gengrid(r) = [vcat(map(x -> x:x:9x, r)...); r[end] * 10]
    p1 = plot(
        x = L1,
        y = L2,
        Geom.line,
        layer(x = mL1, y = mL2, label = label, Geom.point, Geom.label),
        Guide.title("L-curve between λ = 0.01 and 2"),
        Guide.xlabel("Residual norm ||𝐀*𝕟-𝕣||<sub>2</sub>"),
        Guide.ylabel("Solution norm ||𝐈*(𝕟-𝕟ᵢ)||<sub>2</sub>"),
        Guide.xticks(ticks = log10.(gengrid([1]))),
        Guide.yticks(ticks = log10.([gengrid([1, 10]); 200])),
        Theme(
            plot_padding = [0Gadfly.mm, 25Gadfly.mm, 5Gadfly.mm, 0Gadfly.mm],
            default_color = "black",
        ),
        Scale.x_log10(labels = x -> x in ylabels ? @sprintf("%.1f", exp10(x)) : ""),
        Scale.y_log10(labels = x -> x in xlabels ? @sprintf("%i", exp10(x)) : ""),
        Coord.cartesian(
            xmin = log10(1),
            xmax = log10(10),
            ymin = log10(1),
            ymax = log10(200),
        ),
    )

    p2 = plot(
        dfxx,
        x = :Dp,
        y = :S,
        color = :Dist,
        Geom.step,
        Guide.xlabel("Particle diameter (nm)"),
        Guide.ylabel("dN/dlnD (cm-3)"),
        Guide.title("Tikhonov: 𝕟ⁱⁿᵛ=(𝐀ᵀ𝐀+λ²𝐈)⁻¹(𝐀ᵀ𝕣+λ²𝐒⁻¹𝕣)"),
        Guide.xticks(
            ticks = log10.([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500]),
        ),
        Guide.colorkey(; title = ""),
        Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
        Scale.color_discrete_manual(colors...),
        Theme(plot_padding = [0Gadfly.mm, 0Gadfly.mm, 5Gadfly.mm, 0Gadfly.mm]),
        Coord.cartesian(xmin = log10(10), xmax = log10(500)),
    )
    set_default_plot_size(24Gadfly.cm, 8Gadfly.cm)
    p = hstack(p1, p2)
end

lcurvefig()
