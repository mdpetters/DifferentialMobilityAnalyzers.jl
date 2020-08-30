using Distributions, Random, Interact, JLD2

function poisson_noise(Qcpc, 𝕣)

    tscan = 120                # SMPS scan time [s] 
    t = tscan ./ length(𝕣.Dp)    # time in each bin

    c = 𝕣.N * Qcpc * t # number of counts in each bin
    R = Float64[]   # Construct a noisy response function 𝕣+ϵ, where ϵ is counting uncertainty
    for i in c
        f = rand(Poisson(i), 1)   # draw Poisson random number with mean = counts
        push!(R, f[1] / (Qcpc * t))   # convert back to concentration
    end
    R1 = R ./ 𝕣.ΔlnD

    𝕣 = SizeDistribution([], 𝕣.De, 𝕣.Dp, 𝕣.ΔlnD, R ./ 𝕣.ΔlnD, R, :respnse)
end

push!(Base.DL_LOAD_PATH, ".")
function Ωccall1(z, zˢ, β)               # Wrapper for ccall
    return ccall(
        (:omega_, "extlib.so"),
        Float32,
        (Ref{Float32}, Ref{Float32}, Ref{Float32}),
        z,
        zˢ,
        β,
    )
end
Ωext = (Z, zˢ) -> map(z -> Ωccall1(z, zˢ, β), Z);


function compute_matrices(bins)
    t, p = 295.15, 1e5
    qsa, β = 1.66e-5, 1 / 5
    r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369
    leff = 13.0
    m = 6
    DMAtype = :cylindrical
    polarity = :-

    Λ = DMAconfig(t, p, qsa, qsa / β, r₁, r₂, l, leff, polarity, m, DMAtype)
    bins, z₁, z₂ = bins, vtoz(Λ, 10000), vtoz(Λ, 10)
    δ = setupDMA(Λ, z₁, z₂, bins)

    T(zˢ, k, Λ, δ) = δ.Ω(Λ, δ.Z, zˢ / k) .* δ.Tc(k, δ.Dp) .* δ.Tl(Λ, δ.Dp)
    𝐀 = (hcat(map(zˢ -> Σ(k -> T(zˢ, k, Λ, δ), Λ.m), δ.Z)...))'
    𝐀 = 𝐀[:, :]

    df1 = DataFrame(𝐀 = [𝐀], 𝐀⁻¹ = [inv(𝐀)], δ = [δ], Diff = true, Loss = true, bins = bins)

    T(zˢ, k, Λ, δ) = δ.Ω(Λ, δ.Z, zˢ / k) .* δ.Tc(k, δ.Dp)
    𝐀 = (hcat(map(zˢ -> Σ(k -> T(zˢ, k, Λ, δ), Λ.m), δ.Z)...))'
    𝐀 = 𝐀[:, :]
    df2 =
        DataFrame(𝐀 = [𝐀], 𝐀⁻¹ = [inv(𝐀)], δ = [δ], Diff = true, Loss = false, bins = bins)

    T(zˢ, k, Λ, δ) = Ωext(δ.Z, zˢ / k) .* δ.Tc(k, δ.Dp) .* δ.Tl(Λ, δ.Dp)
    𝐀 = (hcat(map(zˢ -> Σ(k -> T(zˢ, k, Λ, δ), Λ.m), δ.Z)...))'
    𝐀 = 𝐀[:, :]
    df3 =
        DataFrame(𝐀 = [𝐀], 𝐀⁻¹ = [inv(𝐀)], δ = [δ], Diff = false, Loss = true, bins = bins)

    T(zˢ, k, Λ, δ) = Ωext(δ.Z, zˢ / k) .* δ.Tc(k, δ.Dp)
    𝐀 = (hcat(map(zˢ -> Σ(k -> T(zˢ, k, Λ, δ), Λ.m), δ.Z)...))'
    𝐀 = 𝐀[:, :]
    df4 =
        DataFrame(𝐀 = [𝐀], 𝐀⁻¹ = [inv(𝐀)], δ = [δ], Diff = false, Loss = false, bins = bins)

    return [df1; df2; df3; df4]
end

function precompute_matrices()
    # Slow
    asd = map([15, 30, 60, 120, 240, 480]) do bins
        compute_matrices(bins)
    end
    xx = vcat(asd...)

    JLD2.jldopen("precomputed_matrices.jld", "w") do file
        write(file, "matrices", xx)
    end
end

function noise_app(Qcpc, N, Bins)
    i =
        (matrices[!, :bins] .== Bins) .&
        (matrices[!, :Diff] .== true) .&
        (matrices[!, :Loss] .== true)
    δ = matrices[i, :δ][1]
    𝐀 = δ.𝐀
    𝐀i = inv(δ.𝐀)

    𝕟 = DMALognormalDistribution([[0.4 * N, 30, 1.2], [0.5 * N, 110, 1.7]], δ)
    lpm = 16.666666 # 16.66 cm³ s⁻¹ = 1 L min⁻¹
    𝕣 = Qcpc <= 1.0 ? poisson_noise(Qcpc * lpm, 𝐀 * 𝕟) : 𝐀 * 𝕟
    𝕟ⁱⁿᵛ = 𝐀i * 𝕣

    df1 = DataFrame(Dp = 𝕟.Dp, S = 𝕟.S, Dist = ["𝕟 (true)" for i = 1:length(𝕟.Dp)])
    df2 = DataFrame(
        Dp = 𝕟ⁱⁿᵛ.Dp,
        S = 𝕟ⁱⁿᵛ.S,
        Dist = ["𝕟ⁱⁿᵛ = 𝐀<sup>-1</sup> 𝕣" for i = 1:length(𝕣.Dp)],
    )
    df = [df1; df2]

    df3 = DataFrame(Dp = 𝕣.Dp, S = 𝕣.N, Dist = ["𝕣" for i = 1:length(𝕣.Dp)])

    xlabels = log10.([10, 20, 50, 100, 200, 500])
    colors = ["black", "darkred"]

    #set_default_plot_size(25Gadfly.cm, 8Gadfly.cm)
    p1 = plot(
        df3,
        x = :Dp,
        y = :S,
        color = :Dist,
        Geom.step,
        Guide.xlabel("Particle diameter (nm)"),
        Guide.ylabel("Concentration (cm-3)"),
        Guide.xticks(
            ticks = log10.([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500]),
        ),
        Guide.colorkey(; title = "Raw \n Response"),
        Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
        Scale.color_discrete_manual(colors...),
        Coord.cartesian(xmin = log10(10), xmax = log10(500)),
    )

    p2 = plot(
        df,
        x = :Dp,
        y = :S,
        color = :Dist,
        Geom.step,
        Guide.xlabel("Particle diameter (nm)"),
        Guide.ylabel("dN/dlnD (cm-3)"),
        Guide.xticks(
            ticks = log10.([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500]),
        ),
        Guide.colorkey(; title = "Inversion"),
        Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
        Scale.color_discrete_manual(colors...),
        Coord.cartesian(xmin = log10(10), xmax = log10(500)),
    )
    set_default_plot_size(24Gadfly.cm, 8Gadfly.cm)
    return hstack(p1, p2)
end

Qcpc = togglebuttons(
    OrderedDict("noise free" => 10.0, "0.05 Lpm" => 0.05, "0.3 Lpm" => 0.3, "1 Lpm" => 1.0),
    value = 0.3,
    label = "CPC Flow",
)

Ntot = togglebuttons(
    OrderedDict("300 cm-3" => 300.0, "3000 cm-3" => 3000.0, "30,000 cm-3" => 30000.0),
    value = 3000.0,
    label = "Total Number Concentration",
)

Bins = togglebuttons(
    OrderedDict(
        "15" => 15,
        "30" => 30,
        "60" => 60,
        "120" => 120,
        "240" => 240,
        "480" => 480,
    ),
    value = 120,
    label = "Bins",
)

Interact.display(hbox(Ntot, Qcpc))
Interact.display(hbox(Bins))

map(noise_app, Qcpc, Ntot, Bins)
