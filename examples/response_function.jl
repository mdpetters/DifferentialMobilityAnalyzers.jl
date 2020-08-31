using DifferentialMobilityAnalyzers
using Gadfly
using DataFrames
using Colors
using Printf

function getresponse(𝕣, 𝕟)
    df1 = DataFrame(Dp = 𝕟.Dp, S = 𝕟.S, Dist = ["𝕟" for i = 1:length(𝕟.Dp)])
    df2 = DataFrame(Dp = 𝕣.Dp, S = 𝕣.S, Dist = ["𝕣" for i = 1:length(𝕣.Dp)])
    df = [df1; df2]

    xlabels = log10.([10, 20, 50, 100, 200, 500])
    colors = ["black", "darkred"]

    set_default_plot_size(16cm, 9cm)
    return plot(
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
        Guide.colorkey(; title = "Distribution"),
        Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
        Scale.color_discrete_manual(colors...),
        Coord.cartesian(xmin = log10(10), xmax = log10(500)),
    )
end

qsa,qsh = 1.66e-5, 8.33e-5           
t,p = 295.15, 1e5                    
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369   
Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,0.0,:-,6,:cylindrical)  
z₁,z₂ = vtoz(Λ,10000), vtoz(Λ,10)    
δ  = setupDMA(Λ, z₁, z₂, 60); 
𝕟 = DMALognormalDistribution([[400, 30, 1.2],[500, 110, 1.7]], δ)
𝕣 = δ.𝐀 * 𝕟
getresponse(𝕣, 𝕟)