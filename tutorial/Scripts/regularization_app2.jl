using Distributions, Random, Interact, JLD2


function regularization_app2(λ)
    Random.seed!(703)
    bins = 120
    i =
        (matrices[!, :bins] .== bins) .&
        (matrices[!, :Diff] .== true) .&
        (matrices[!, :Loss] .== true)
    δ = matrices[i, :δ][1]
    𝐀 = δ.𝐀
    𝐒 = δ.𝐒
    𝐒i = inv(δ.𝐒)
    𝐀i = inv(δ.𝐀)
    N = 500.0
    𝕟 = DMALognormalDistribution([[0.4 * N, 30, 1.2], [0.5 * N, 110, 1.7]], δ)
    lpm = 16.666666 # 16.66 cm³ s⁻¹ = 1 L min⁻¹
    𝕣 = poisson_noise(1 * lpm, 𝐀 * 𝕟)
    𝕟ⁱⁿᵛ = 𝐀i * 𝕣

    df1 = DataFrame(Dp = 𝕟.Dp, S = 𝕟.S, Dist = ["𝕟 (true)" for i = 1:length(𝕟.Dp)])
    df2 = DataFrame(Dp = 𝕟ⁱⁿᵛ.Dp, S = 𝕟ⁱⁿᵛ.S, Dist = ["𝕟ⁱⁿᵛ" for i = 1:length(𝕣.Dp)])
    df = [df1; df2]

    𝕟ⁱⁿᵛ = 𝐒i * 𝕣
    df1 = DataFrame(Dp = 𝕟.Dp, S = 𝕟.S, Dist = ["𝕟 (true)" for i = 1:length(𝕟.Dp)])
    df2 = DataFrame(
        Dp = 𝕟ⁱⁿᵛ.Dp,
        S = [0.0 for i = 1:length(𝕟ⁱⁿᵛ.Dp)],
        Dist = ["𝕟ⁱⁿᵛ=0" for i = 1:length(𝕟ⁱⁿᵛ.Dp)],
    )
    dfx = [df1; df2]

    𝐈 = Matrix{Float64}(I, bins, bins)
    #𝕟ⁱⁿᵛ = inv(𝐀'*𝐀 + λ^2.0*𝐈) * (𝐀'𝕣  + λ^2.0 * (inv(𝐒)*𝕣))
    if λ > 0
        #𝕟ⁱⁿᵛ = (𝐀'𝐀 + λ^2𝐈)^(-1) * (𝐀'𝕣  + λ^2 * 𝐒^(-1)*𝕣)
        𝕟ⁱⁿᵛ = (𝐀'𝐀 + λ^2𝐈)^(-1) * (𝐀'𝕣)
    else
        𝕟ⁱⁿᵛ = (𝐀)^(-1) * 𝕣
    end

    df1 = DataFrame(Dp = 𝕟.Dp, S = 𝕟.S, Dist = ["𝕟 (true)" for i = 1:length(𝕟.Dp)])
    df2 = DataFrame(Dp = 𝕟ⁱⁿᵛ.Dp, S = 𝕟ⁱⁿᵛ.S, Dist = ["𝕟ⁱⁿᵛ" for i = 1:length(𝕟ⁱⁿᵛ.Dp)])
    dfxx = [df1; df2]

    df3 = DataFrame(Dp = 𝕣.Dp, S = 𝕣.N, Dist = ["𝕣" for i = 1:length(𝕣.Dp)])

    xlabels = log10.([10, 20, 50, 100, 200, 500])
    colors = ["black", "darkred"]

    set_default_plot_size(28Gadfly.cm, 8Gadfly.cm)
    p1 = plot(
        df,
        x = :Dp,
        y = :S,
        color = :Dist,
        Geom.step,
        Guide.xlabel("Particle diameter (nm)"),
        Guide.ylabel("Concentration (cm-3)"),
        Guide.xticks(
            ticks = log10.([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500]),
        ),
        Guide.title("Least-squares: 𝕟ⁱⁿᵛ = 𝐀<sup>-1</sup> 𝕣"),
        Guide.colorkey(; title = ""),
        Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
        Scale.color_discrete_manual(colors...),
        Coord.cartesian(xmin = log10(10), xmax = log10(500)),
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
        Coord.cartesian(xmin = log10(10), xmax = log10(500)),
    )


    p3 = plot(
        dfx,
        x = :Dp,
        y = :S,
        color = :Dist,
        Geom.step,
        Guide.xlabel("Particle diameter (nm)"),
        Guide.ylabel("dN/dlnD (cm-3)"),
        Guide.xticks(
            ticks = log10.([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500]),
        ),
        Guide.title("Initial Guess: 𝕟ⁱⁿᵛ = 0"),
        Guide.colorkey(; title = ""),
        Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
        Scale.color_discrete_manual(colors...),
        Coord.cartesian(xmin = log10(10), xmax = log10(500)),
    )

    return hstack(p1, p2, p3)
end

myλ = slider([0.0; 0.01:0.01:0.09; 0.1:0.1:0.9; 1.0:1.0:10.0], value = 0.01, label = "λ")

Interact.display(myλ)

map(regularization_app2, myλ)
