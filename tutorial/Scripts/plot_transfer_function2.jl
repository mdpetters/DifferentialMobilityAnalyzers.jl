using Gadfly, DifferentialMobilityAnalyzers, DataFrames, Colors

function getΩ2()
    t, p = 295.15, 1e5
    qsa, β = 1.66e-5, 1 / 5
    r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369
    leff = 13.0
    m = 6
    DMAtype = :cylindrical
    polarity = :-

    Λ = DMAconfig(t, p, qsa, qsa / β, r₁, r₂, l, leff, polarity, m, DMAtype)
    bins, z₁, z₂ = 5, vtoz(Λ, 10000), vtoz(Λ, 10)
    δ = setupDMA(Λ, z₁, z₂, bins)

    tscan = 60                 # 60 second SMPS scan
    tc = 2.0                   # 2 second integration time per channel
    vmin, vmax = 10, 10000       # Scan from 10V to 10kV
    bins = round(Int, tscan / tc) # Number of size bins
    Ve = reverse(exp10.(range(log10(vmin), stop = log10(vmax), length = bins + 1)))  # Voltage bin-edges, down scan 
    Vp = sqrt.(Ve[2:end] .* Ve[1:end-1])  # Voltage midpoints
    t = 1:tc:tscan                      # Time array
    Z = exp10.(range(log10(1e-9), stop = log10(1e-6), length = 1000))


    Zˢ = vtoz(Λ, Ve)     # Vector of mobility centroids for bin edges
    i = 10               # Bin number to obtain average transfer function
    nint = 20             # Number of points for numerical integration


    # Numerical integration over all transfer functions in the bin
    Vex = reverse(exp10.(range(log10(Ve[i+1]), stop = log10(Ve[i]), length = nint)))
    Ωav = mapreduce(zˢ -> δ.Ω(Λ, Z, zˢ), +, vtoz(Λ, Vex)) / nint


    #plot(Z, ), ylabel = "Ω (-)", ylim = (0,1), xlabel = "Z (m² V⁻¹ s⁻¹)", color = :black,
    #left_margin = 33px, label = "Ω @ tᵢ", xaxis = :log10, xlim = (10e-9, 50e-9))
    #plot!(Z, δ.Ω(Λ,Z,Zˢ[i+1]), label = "Ω @ tᵢ+tc", color = :black, ls = :dashdot)
    #plot!(Z,Ωav, color = RGBA(0.8,0,0,1), label = "Ωav",fmt=:svg)


    zˢ = Z[argmax(Ωav)]

    df1 = DataFrame(
        x = Z ./ zˢ,
        Ω = δ.Ω(Λ, Z, Zˢ[i]),
        Dp = ["Ω (-) @ tᵢ" for i = 1:length(Z)],
    )
    df2 = DataFrame(
        x = Z ./ zˢ,
        Ω = δ.Ω(Λ, Z, Zˢ[i+1]),
        Dp = ["Ω @ tᵢ+tc" for i = 1:length(Z)],
    )
    df3 = DataFrame(x = Z ./ zˢ, Ω = Ωav, Dp = ["Ωav" for i = 1:length(Z)])
    df = [df1; df2; df3]
    set_default_plot_size(16cm, 9cm)

    p = plot(
        df,
        x = :x,
        y = :Ω,
        color = :Dp,
        Geom.line(preserve_order = true),
        Guide.xlabel("z/z <sup>s</sup> (-)"),
        Guide.ylabel("Ω (-)"),
        Guide.colorkey(; title = ""),
        Guide.xticks(ticks = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]),
        Guide.yticks(ticks = 0:0.2:1),
        Scale.color_discrete_manual(["darkred", "steelblue3", "black"]...),
        Coord.cartesian(xmin = 0.6, xmax = 1.4, ymin = 0, ymax = 1),
    )
end

getΩ2()
