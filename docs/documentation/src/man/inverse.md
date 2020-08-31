# Inverse Models

## Size Distribution Inversion

### Theory
The size distribution can be found from the inverse of the response function using [Tikhonov regularization](https://en.wikipedia.org/wiki/Tikhonov_regularization).

``L_1`` and ``L_2`` are the Euclidian norms of the solution and the intial guess. 

``L_1 = \left\lVert \bf{A}\rm{𝕟} - \rm{𝕣}\right\rVert_2``

``L_2 = \left\lVert\bf{L}(\rm{𝕟} - \rm{𝕟_i})\right\rVert_2``

where ``\rm{𝕟}`` is the true size distribution, ``\bf{A}`` is the convolution matrix, 𝕣 is the measured response vector, ``\bf{L}`` is a weights matrix, and ``\rm{𝕟_i}`` is an initial guess. To solve for the inverse, a balance is sought between the ``L_1`` and ``L_2`` norms.

``𝕟_{inv} = \arg \min\{L_1^2 + \lambda^2 L_2^2\}``

where ``\lambda`` is the regularization parameter. Taking ``\bf{L} = \bf{I}`` (weight matrix equals the identity [Matrix 𝐈](@ref)) and ``\rm{𝕟_i} = \bf{S}^{-1}\rm{𝕣}`` as initial guess ([Matrix 𝐒](@ref)), the regularized inverse is computed via:

``𝕟_{inv} = (\bf{A}^{\rm{T}}\bf{A} + \lambda^{\rm{2}} \bf{I})^{\rm{-1}}(\bf{A}^{\rm{T}} \rm{𝕣} + \lambda^2\bf{S}^{-1} \rm{𝕣})``

The regularization parameter ``\lambda`` "interpolates" between the noisy least-square inverse (``\lambda = 0``) and the initial guess (``\lambda >> 1``). The optimal regularization parameter is found using the L-curve method. The optimal ``\lambda_{opt}`` is found using the L-curve method. The L-curve is defined as a plot of ``\log(L_1)`` vs.  ``\log(L_2)``, where ``L_1`` and ``L_2`` are obtained for a series of discrete ``\lambda`` values. The ``\lambda_{opt}`` is found at the corner of the L-curve, which is mathematically defined as the point where where the curvature of the L-curve is maximum. Here, the corner of the L-curve is found using the iterative algorithm described in Talukdar and Swihart (2003). The curvature is calculated using Eq. (14) in Hansen (2000), which requires the first and second derivatives of ``d\ln(L_i)^2/d\lambda``. These derivatives of ``d\ln(L_i)^2/d\lambda`` are estimated numerically.

See session 3 of the [Tutorial](@ref) for an interactive explanation on how Tikhonov regularization works.

### Example
The L-curve search method is implemented in the [rinv](@ref) function. 

```julia
𝕟ⁱⁿᵛ = rinv(𝕣.N, δ, λ₁=0.1, λ₂=1.0);
```

The function takes the response vector to be inverted, the DMA grid to read the matrices, and an upper and lower search bound for the optimal regularization parameter. 

```@example
using DataFrames, Gadfly, CSV, Printf, DifferentialMobilityAnalyzers #hide
#hide
df = CSV.read("example_data.csv", DataFrame) 

t, p, lpm = 293.15, 940e2, 1.666e-5 
r₁, r₂, l = 9.37e-3,1.961e-2,0.44369     
Λ = DMAconfig(t,p,1lpm,4lpm,r₁,r₂,l,0.0,:+,6,:cylindrical)  
δ = setupDMA(Λ, vtoz(Λ,10000), vtoz(Λ,10), 120) 

𝕣 = (df,:Dp,:Rcn,δ) |> interpolateDataFrameOntoδ 
𝕟ⁱⁿᵛ = rinv(𝕣.N, δ, λ₁=0.1, λ₂=1.0)

df = DataFrame(Dp = 𝕟ⁱⁿᵛ.Dp, S = 𝕟ⁱⁿᵛ.S, Dist = ["𝕟ⁱⁿᵛ" for i = 1:length(𝕟ⁱⁿᵛ.Dp)]) #hide
dfr = DataFrame(Dp = 𝕣.Dp, S = 𝕣.N, Dist = ["𝕣" for i = 1:length(𝕣.Dp)])#hide
#hide
xlabels = log10.([10, 100, 500])#hide
colors = ["darkred", "steelblue3", "black"]#hide
p1 = plot(#hide
    dfr,#hide
    x = :Dp,#hide
    y = :S,#hide
    color = :Dist,#hide
    Geom.step,#hide
    Guide.xlabel("Particle diameter (nm)"),#hide
    Guide.ylabel("N (cm-3)"),#hide
    Guide.title("Raw Response Function"),#hide
    Guide.colorkey(; title = ""),#hide
    Guide.xticks(#hide
        ticks = log10.([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500]),#hide
    ),#hide
    Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),#hide
    Scale.color_discrete_manual("black"),#hide
    Coord.cartesian(xmin = log10(10), xmax = log10(500), ymin = 0),#hide
)#hide
#hide
p2 = plot(#hide
    df,#hide
    x = :Dp,#hide
    y = :S,#hide
    color = :Dist,#hide
    Geom.step,#hide
    Guide.xlabel("Particle diameter (nm)"),#hide
    Guide.ylabel("dN/dlnD (cm-3)"),#hide
    Guide.title("Inverted Size Distribution"),#hide
    Guide.xticks(#hide
        ticks = log10.([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500]),#hide
    ),#hide
    Guide.colorkey(; title = ""),#hide
    Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),#hide
    Scale.color_discrete_manual(colors...),#hide
    Coord.cartesian(xmin = log10(10), xmax = log10(500), ymin = 0),#hide
)#hide
set_default_plot_size(20cm, 8cm)#hide
p = (hstack(p1, p2))#hide
```

The file ```example_data.csv``` is a comma delimited text files that contains diameter and a resonse vector. The function [interpolateDataFrameOntoδ](@ref) 
takes the columns :Dp and :Rcn from the DataFrame df and creates a response distribution of the type SizeDistribution that matches the DMA grid δ. Finally, 𝕟ⁱⁿᵛ is computed using the [rinv](@ref) function. 

The method was compared to the inverted size distribution given by the
manufacturer Aerosol Instrument Manager (AIM) software, assuming no diffusion
correction in either method. Figure 3 in the [manuscript](https://www.tandfonline.com/doi/full/10.1080/02786826.2018.1530724) shows excellent correspondence between this package and the commercial software.

## CCN Inversion

Size resolved CCN measurements have a long history (e.g. Cruz and Pandis, 1997, Snider et al., 2006). Particles are typically dried, charge neutralized, and passed through the DMA.At the exit, the flow is split between a CPC that measures all particles and a CCN counter that measures particles that form cloud droplets at a water supersaturation set by the instrument. In this configuration the DMA voltage can either be stepped or continuously scanned. The ratio of CCN-to-CPC response function is used to determine the fraction of particles that activate at a given diameter. The diameter where 50% of the particles activate is taken to be the activation diameter. The activation of particles into cloud droplets is proportional to the volume of solute in the particle. Therefore, larger multiply charged particles activate first. This leads to a bias in the inferred D50 diameter if the activated fraction is computed from the ratio of the raw response functions (Petters et al., 2007, 2009).

The CCN transmission function is modeled using a cumulative Gauss integral

``T_{af} = \frac{1}{2}\left[1 + \mathrm{erf}\left(\frac{d-\mu}{\sigma}\right) \right]``

This function is applied to the mobility size distribution. Then the response function is computed. Empirically, activated fraction can be computed using from the ratio of size distributions and response functions. 

To predict the ratio of 𝕣ᶜᶜⁿ/𝕣ᶜⁿ in terms of the activated fraction model Taf:

``𝕟^{\mathrm{cn}} = (\bf{A}^{\rm{T}}\bf{A} + \lambda^{\rm{2}} \bf{I})^{\rm{-1}}(\bf{A}^{\rm{T}} \rm{𝕣^{\mathrm{cn}}} + \lambda^2 \bf{S}^{-1} \rm{𝕣^{\mathrm{cn}}})``

``\left( \frac{𝕣^{\mathrm{ccn}}}{𝕣^{\mathrm{cn}}} \right)_{\mathrm{model}} = \frac{\mathbf{A}[T_{af}\; .*\; 𝕟^{\mathrm{cn}}]}{\mathbf{A}𝕟^{\mathrm{cn}}}``

This mirrors the approach in Petters et al. (2007). The model + fit can be solved using this package as follows.

```julia
𝕟ᶜⁿ = (𝐀'𝐀 + λ^2𝐈)^(-1) * (𝐀'𝕣ᶜⁿ  + λ^2 * 𝐒^(-1)*𝕣ᶜⁿ)
model(x,p) = (δ.𝐀*(𝕟ᶜⁿ.N.*Taf(𝕟ᶜⁿ, p[1], p[2])))./(δ.𝐀*𝕟ᶜⁿ.N)
fit = curve_fit(model, δ.Dp, 𝕣ᶜᶜⁿ.N./𝕣ᶜⁿ.N, [65.0, 3.0])
```
!!! note
    The example below includes a hidden threholding step. (See complete example in the examples folder) Thresholding is performed for bins with low counts, which results in NaN and Inf in 𝕒𝕗. This is necessary because at low concentrations a bin may have zero or very low concentration, resulting in unrealistic or NaN/InF values in the activated fraction array. Thus if counts/concentration is below the threshold, the activated fraction is forced to one for large diameters and zero for small diameters. The noise threshold may need to be adjusted for different datasets.

The code loads the CCN (𝕣ᶜᶜⁿ) and CN (𝕣ᶜⁿ) response functions from file. The activated fraction is computed through the ratio of the distributions using one of the SizeDistribution [Number Operators](@ref).

```julia
𝕒𝕗 = 𝕣ᶜᶜⁿ/𝕣ᶜⁿ
```

And here is the complete example showing how to fit the 𝕒𝕗 response curve.

```@example
using DataFrames, Gadfly, CSV, DifferentialMobilityAnalyzers, Printf, LsqFit, SpecialFunctions #hide
t, p, lpm = 293.15, 940e2, 1.666e-5     
r₁, r₂, l = 9.37e-3,1.961e-2,0.44369    
Λ = DMAconfig(t,p,1lpm,4lpm,r₁,r₂,l,0.0,:+,6,:cylindrical)  
δ = setupDMA(Λ, vtoz(Λ,10000), vtoz(Λ,10), 120)

# Load a simple comma delimited text file - file contains :Dp, :Rcn, :Rccn
df = CSV.read("example_data.csv")
𝕣ᶜⁿ = (df,:Dp,:Rcn,δ) |> interpolateDataFrameOntoδ        # CN response distribution
𝕣ᶜᶜⁿ = (df,:Dp,:Rccn,δ) |> interpolateDataFrameOntoδ;     # CCN response distribution
#hide
function threshold!(𝕟::SizeDistribution, c::Float64, n1::Float64, n2::Float64)#hide
    N = 𝕟.N  #hide
    S = 𝕟.S#hide
    S[(N .<= c) .& (𝕟.Dp .> 150)] .= n2#hide
    N[(N .<= c) .& (𝕟.Dp .> 150)] .= n2#hide
    S[(N .<= c) .& (𝕟.Dp .< 150)] .= n1#hide
    N[(N .<= c) .& (𝕟.Dp .< 150)] .= n1#hide
    𝕟.N = N#hide
 end#hide
 #hide
threshold!(𝕣ᶜⁿ, 0.1, 0.1, 0.1)#hide
threshold!(𝕣ᶜᶜⁿ, 0.1, 0.0, 0.1)#hide
#hide
𝕒𝕗 = 𝕣ᶜᶜⁿ/𝕣ᶜⁿ
Taf(𝕟,μ,σ) = @. 0.5 * (1.0 + erf((𝕟.Dp - μ)./(sqrt(2.0σ))))
𝐈, 𝐒, 𝐀, λ =  δ.𝐈, δ.𝐒, δ.𝐀, 0.5
𝕟ᶜⁿ = (𝐀'𝐀 + λ^2𝐈)^(-1) * (𝐀'𝕣ᶜⁿ  + λ^2 * 𝐒^(-1)*𝕣ᶜⁿ)
model(x,p) = (𝐀 * (𝕟ᶜⁿ.N .* Taf(𝕟ᶜⁿ, p[1], p[2])))./( 𝐀 * 𝕟ᶜⁿ.N)
fit = curve_fit(model, 𝕒𝕗.Dp, 𝕒𝕗.N, [65.0, 3.0])
Ax = fit.param
afmodel = model(δ.Dp, Ax)
#hide
df1 = DataFrame(Dp = 𝕣ᶜⁿ.Dp, S = 𝕣ᶜⁿ.S, Dist = ["𝕣ᶜⁿ" for i = 1:length(𝕣ᶜⁿ.Dp)])#hide
df2 = DataFrame(Dp = 𝕣ᶜᶜⁿ.Dp, S = 𝕣ᶜᶜⁿ.S, Dist = ["𝕣ᶜᶜⁿ" for i = 1:length(𝕣ᶜᶜⁿ.Dp)])#hide
df = [df1; df2]#hide
#hide
dfr1 = DataFrame(Dp = 𝕒𝕗.Dp, S = 𝕒𝕗.N, Dist = ["𝕒𝕗 (data)" for i = 1:length(𝕒𝕗.Dp)])#hide
dfr2 = DataFrame(Dp = 𝕒𝕗.Dp, S = afmodel, Dist = ["𝕒𝕗 (model)" for i = 1:length(𝕒𝕗.Dp)]) #hide
dfr = [dfr1; dfr2]#hide
#hide
xlabels = log10.([10, 100, 500])#hide
colors = ["darkred", "steelblue3", "black"]#hide
p1 = plot(#hide
    dfr,#hide
    x = :Dp,#hide
    y = :S,#hide
    color = :Dist,#hide
    Geom.step,#hide
    Guide.xlabel("Apparent +1 Mobility Diameter (nm)"),#hide
    Guide.ylabel("Fraction (-)"),#hide
    Guide.title("Activated Fraction"),#hide
    Guide.colorkey(; title = ""),#hide
    Guide.xticks(#hide
        ticks = log10.([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500]),#hide
    ),#hide
    Guide.yticks(ticks = collect(0:0.2:1.2)),#hide
    Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),#hide
    Scale.color_discrete_manual(["black", "darkgoldenrod"]...),#hide
    Coord.cartesian(xmin = log10(10), xmax = log10(500), ymin = 0, ymax = 1.2),#hide
)#hide
#hide
p2 = plot(#hide
    df,#hide
    x = :Dp,#hide
    y = :S,#hide
    color = :Dist,#hide
    Geom.step,#hide
    Guide.xlabel("Apparent +1 Mobility Diameter (nm)"),#hide
    Guide.ylabel("dN/dlnD (cm-3)"),#hide
    Guide.title("Raw response function"),#hide
    Guide.xticks(#hide
        ticks = log10.([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500]),#hide
    ),#hide
    Guide.colorkey(; title = ""),#hide
    Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),#hide
    Scale.color_discrete_manual(colors...),#hide
    Coord.cartesian(xmin = log10(10), xmax = log10(500), ymin = 0),#hide
)#hide
set_default_plot_size(20cm, 8cm)#hide
p = (hstack(p2, p1))#hide
```