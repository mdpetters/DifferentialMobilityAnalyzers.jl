# Forward Models

## Single DMA

### Constant Voltage
The DMA is operated in classifier mode and set to a single voltage. The mobility classified particles are then passed to one or more detectors. Transmission can be modeled by the convolution of the particle [Charging Probability](@ref), DMA [Transfer Function](@ref), and [Transmission Loss](@ref) through the DMA. These are built into the [DMA Grid](@ref). Therefore a net transmission function can be written:

!!! info
    The net transmission function T describes transmission through the DMA, including the charge filter, the DMA transfer function, and the DMA transmission efficiency. 
    ```julia
    T(zˢ,k,Λ,δ) = δ.Ω(Λ,δ.Z,zˢ/k).*δ.Tc(k,δ.Dp).*δ.Tl(Λ,δ.Dp)
    ```
    where Λ is the DMA configuration, δ the DMA grid, zˢ is the z-star selected by the DMA, and k is the number of charges. 	

The following example shows ho to use T(zˢ,k,Λ,δ) to predict the mobility distribution exiting the DMA. In the code fragment
- zˢ is the z-star selected by the DMA
- 𝕟ᶜⁿ is an assumed known bimodal lognormal distribution computed on the DMA size grid
- ℕ is an array of transmitted mobility distributions carrying k charges
- 𝕄 is an array of transmitted apparent mobility distributions carrying k charges
- 𝕟ₜ, 𝕞ₜ are the superposition of these distributions


```@example
using Gadfly, NumericIO, Colors, LinearAlgebra, Printf, DataFrames, Cairo, Fontconfig # hide
using DifferentialMobilityAnalyzers

# Create a DMA config
qsa,qsh = 1.66e-5, 8.33e-5                       # Qsample [m3 s-1], Qsheath [m3 s-1]
t,p = 295.15, 1e5                                # Temperature [K], Pressure [Pa]
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369               # DMA geometry [m]
Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,0.0,:-,6,:cylindrical)  

# Create a DMA grid
z₁,z₂ = vtoz(Λ,10000), vtoz(Λ,10)    # bins, upper, lower mobility limit
δ = setupDMA(Λ, z₁, z₂, 512); 

# Compute the transmission through the DMA
T(zˢ,k,Λ,δ) = δ.Ω(Λ,δ.Z,zˢ/k).*δ.Tc(k,δ.Dp).*δ.Tl(Λ,δ.Dp)                  
zˢ = dtoz(Λ, 100*1e-9)                                                    
𝕟ᶜⁿ = DMALognormalDistribution([[900., 40., 1.5], [500., 180., 1.4]], δ)   
ℕ = map(k -> T(zˢ,k,Λ,δ)*𝕟ᶜⁿ,1:3)                                   
𝕄 = map(k -> (ztod(Λ,1,zˢ)/ztod(Λ,k,zˢ))⋅(T(zˢ,k,Λ,δ)*𝕟ᶜⁿ),1:3)    
𝕟ₜ, 𝕞ₜ = sum(ℕ), sum(𝕄)
# hide
set_default_plot_size(22cm, 6cm) # hide
# hide
xlabels = log10.([10, 50, 100, 500]) # hide
p1 = plot( # hide
    x = 𝕟ᶜⁿ.Dp, # hide
    y = 𝕟ᶜⁿ.S, # hide
    Geom.step, # hide
    color = ["𝕟ᶜⁿ" for i in 𝕟ᶜⁿ.Dp], # hide
    Guide.xlabel("Particle diameter (nm)", orientation = :horizontal), # hide
    Guide.ylabel("dN/dlnD (cm-3)", orientation = :vertical), # hide
    Guide.xticks( # hide
        ticks = log10.([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600]), # hide
    ), # hide
    Guide.colorkey(; title = ""),# hide
    Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),# hide
    Scale.color_discrete_manual("black"),# hide
    Coord.cartesian(xmin = log10(10), xmax = log10(600)),# hide
    Theme(plot_padding = [0mm, 8mm, 0mm, 0mm]),# hide
) # hide
# hide
df1 = let# hide
    xx = map(1:3) do i# hide
        df = DataFrame(x = 𝕄[i].Dp, y = 𝕄[i].S, c = ["𝕄[$i]" for j in 𝕄[i].Dp])# hide
    end# hide
    vcat(xx...)# hide
end# hide
# hide
df2 = DataFrame(x = 𝕞ₜ.Dp, y = 𝕞ₜ.S, c = ["𝕞ₜ" for j in 𝕞ₜ.Dp])# hide
# hide
df = [df2; df1]# hide
# hide
colors = ["black", "darkred", "steelblue3", "darkgoldenrod"]# hide
# hide
p2 = plot(# hide
    df,# hide
    x = :x,# hide
    y = :y,# hide
    color = :c,# hide
    Geom.line,# hide
    Guide.xlabel("Apparent +1 Mobility Diameter (nm)", orientation = :horizontal),# hide
    Guide.ylabel(""),# hide
    Guide.xticks(ticks = [80, 100, 120, 140]),# hide
    Guide.colorkey(; title = ""),# hide
    Scale.color_discrete_manual(colors...),# hide
    Coord.cartesian(xmin = 70, xmax = 150),# hide
    Theme(plot_padding = [5mm, 10mm, 0mm, 0mm]),# hide
)# hide
# hide
df1 = let# hide
    xx = map(1:3) do i# hide
        df = DataFrame(x = ℕ[i].Dp, y = ℕ[i].S, c = ["ℕ[$i]" for j in ℕ[i].Dp])# hide
    end# hide
    vcat(xx...)# hide
end# hide
# hide
df2 = DataFrame(x = 𝕟ₜ.Dp, y = 𝕟ₜ.S, c = ["𝕟ₜ" for j in 𝕟ₜ.Dp])# hide
# hide
df = [df2; df1]# hide
# hide
p3 = plot(# hide
    df,# hide
    x = :x,# hide
    y = :y,# hide
    color = :c,# hide
    Geom.line,# hide
    Guide.xlabel("Mobility Diameter (nm)", orientation = :horizontal),# hide
    Guide.ylabel(""),# hide
    Guide.xticks(ticks = [100, 150, 200, 250]),# hide
    Guide.colorkey(; title = ""),# hide
    Scale.color_discrete_manual(colors...),# hide
    Coord.cartesian(xmin = 70, xmax = 250),# hide
    Theme(plot_padding = [0mm, 10mm, 0mm, 0mm]),# hide
)# hide
# hide
hstack(p1, p2, p3)# hide
```
**Figure.** Left: assumed bimodal lognormal size distribution. Middle: monodisperse mobility size distribution plotted against the apparent +1 mobility diameter, defined as the apparent setpoint diameter of the DMA. Dashed line is total number concentration. Right: same as middle panel but plotted versus the mobility diameter.

The DMA selects approximately a triangular distribution around mobility centroid ``z^s``. The distribution is symmetric when plotted against the log of mobility, but asymmetric when plotted against the log of the apparent +1 mobility diameter because the Cunningham slip flow correction factor applied in the conversion from mobility to diameter is a strong function of particle size. The majority of selected particles are singly charged, but the contribution of multiply charged particles to the total number is not negligible. The mobility distribution has contributions from particles that are at least twice the diameter of the selected centroid diameter. The relative fractions are determined by the equilibrium charge fraction and the number of particles available at each diameter.

!!! info    
    The ```map(k->f, 1:3)``` construct sequentially applies values from the array [1,2,3] to k and calls the function f. The output is an array of length 3. Since ```T(zˢ,k,Λ,δ)``` produces a vector, 𝕟ᶜⁿ is a size distribution, and vector * size distribution is a size distribution, the output of
    ```
    ℕ = map(k -> T(zˢ,k,Λ,δ)*𝕟ᶜⁿ,1:3) 
    ```
    is an array of size distributions.

For more information see Figure 2 in the [Manuscript](https://www.tandfonline.com/doi/full/10.1080/02786826.2018.1530724), check out Session 2 of the [Tutorial](@ref) and/or Notebook S4 in the [Notebooks](@ref) section.

### Complete Distribution

If the DMA is stepped or scanned, the response function can be computed from the convolution [Matrix 𝐀](@ref) and the known true size distribution. The operation matrix * size distribution is defined in the [Operators](@ref) section. Thus we can conveniently write ```𝕣 = 𝐀 * 𝕟``` to compute the response function.

```@example
using DifferentialMobilityAnalyzers #hide
using Gadfly #hide
using DataFrames #hide
using Colors #hide
using Printf #hide
qsa,qsh = 1.66e-5, 8.33e-5               #hide
t,p = 295.15, 1e5                       #hide
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369         #hide
Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,0.0,:-,6,:cylindrical)  #hide
#hide
z₁,z₂ = vtoz(Λ,10000), vtoz(Λ,10)    # #hide
δ = setupDMA(Λ, z₁, z₂, 60); 
#hide
𝕟 = DMALognormalDistribution([[400, 30, 1.2],[500, 110, 1.7]], δ)
𝕣 = δ.𝐀 * 𝕟
#hide
function getresponse(𝕣, 𝕟)#hide
    df1 = DataFrame(Dp = 𝕟.Dp, S = 𝕟.S, Dist = ["𝕟" for i = 1:length(𝕟.Dp)])#hide
    df2 = DataFrame(Dp = 𝕣.Dp, S = 𝕣.S, Dist = ["𝕣" for i = 1:length(𝕣.Dp)])#hide
    df = [df1; df2]#hide
#hide
    xlabels = log10.([10, 20, 50, 100, 200, 500])#hide
    colors = ["black", "darkred"]#hide
#hide
    set_default_plot_size(16cm, 9cm)#hide
    return plot(#hide
        df,#hide
        x = :Dp,#hide
        y = :S,#hide
        color = :Dist,#hide
        Geom.step,#hide
        Guide.xlabel("Particle diameter (nm)"),#hide
        Guide.ylabel("dN/dlnD (cm-3)"),#hide
        Guide.xticks(#hide
            ticks = log10.([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500]),#hide
        ),#hide
        Guide.colorkey(; title = "Distribution"),#hide
        Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),#hide
        Scale.color_discrete_manual(colors...),#hide
        Coord.cartesian(xmin = log10(10), xmax = log10(500)),#hide
    )#hide
end#hide
getresponse(𝕣, 𝕟) # hide
```

## Tandem DMAs

### Humidified Tandem DMA

#### Single Composition
Dried, charge equilibrated particles are classified in DMA1. The flow is split to measure particle concentration with a condensation particle counter (CPC). The remaining flow is passed through a humidifier. Hygroscopic particles take up water and increase in diameter. The humidified size distribution is measured using the second DMA that is operated in scanning or stepping mode. Passage through a second bipolar charger (charge neutralizer) is optional and rarely used in TDMA experiments.

To model transmission through the tandem DMA we need to 
- setup two DMAs, δ₁ and δ₂ 
- know the input size distribution 
- formulate a transmission model

Setting up DMAs is described in [DMA Configuration](@ref). The input size distribution is assumed based on a lognormal distribution. The constructor function [DMALognormalDistribution](@ref) initializes a lognormal distribution along a DMA grid.

```julia
Ax = [[1300.0, 60.0, 1.4], [2000.0, 200.0, 1.6]]
𝕟ᶜⁿ = DMALognormalDistribution(Ax, δ₁)
```

The transmission model is a combination of operating DMA₁ at [Constant Voltage](@ref) and transmission of a [Complete Distribution](@ref) through DMA₂. 

```julia
# Tandem DMA equations
O(k) = mapfoldl(zs -> (δ₂.Ω(Λ₂, δ₂.Z, zs / k, k) .* δ₂.Tl(Λ₂, δ₂.Z, k))', vcat, δ₂.Z)
T₁(zˢ, k) = δ₁.Ω(Λ₁, δ₁.Z, zˢ / k, k) .* δ₁.Tc(k, δ₁.Dp) .* δ₁.Tl(Λ₁, δ₁.Dp)
DMA₁(𝕟, zˢ, gf) = @_ map((gf ⋅ (T₁(zˢ, _) * 𝕟)), 1:6)
itp(𝕟) = interpolateSizeDistributionOntoδ((𝕟, δ₂))
DMA₂(𝕟, k) = O(k) * 𝕟
```

The function ```T(zˢ, k, Λ, δ)``` is already known. The function ```DMA₁(𝕟, zˢ, gf)``` takes a distribution 𝕟 and mobility zˢ and passes it through DMA Λ₁, δ₁ and applied growth factor gf. The resulting distributions are interpolated into the same grid as DMA2 using
[interpolateSizeDistributionOntoδ](@ref). 
The function ```DMA₂(𝕟, δ)``` takes an input size distribution 𝕟 and passes it through DMA₂. No neutralizer is used. Therefore the convolution O(k) is applied. Note that O(k) corresponds to Eq. (15) in [Petters (2021)](https://amt.copernicus.org/articles/14/7909/2021/amt-14-7909-2021.pdf).  

!!! note
    The dot product of scalar ⋅ SizeDistribution shifts the size distribution in diameter space: [Size Operators](@ref). Check out the [Tutorial](@ref) Session 1 and/or Notebook S3 in the [Notebooks](@ref) section for visualizations.

Here is an abriged example how to compute the grown output distributions from DMA2

```julia
Dd = 100e-9             # Dry diameter
zˢ = dtoz(Λ₁, Dd);      # Mobility of 100 nm particle
gf = 1.6                # Growth factor
ℕ = DMA₁(𝕟ᶜⁿ, zˢ, gf)   # Transmission through DMA1
𝕄 = map(k -> (@> itp(ℕ[k]) DMA₂(k)), 1:3) # Transmission through DMA2
𝕞ᵗ = sum(𝕄)                               # total response
```

𝕄[k] correspond to the +1, +2, +3 partial mobility response functions
that would be measured after DMA2. The total is obtained by the sum of these distributions
Below is the complete example, which produces the response function of the tandem DMA.

```@example
using DifferentialMobilityAnalyzers #hide
using Gadfly #hide
using NumericIO #hide
using Colors #hide
using LinearAlgebra #hide
using Printf #hide
using DataFrames #hide
using Underscores #hide
import Lazy.@>, Lazy.@>>#hide
t, p = 295.15, 1e5                             # Temperature [K], Pressure [Pa]
qsa, β = 1.66e-5, 1 / 5                        # Qsample [m3 s-1], Sample-to-sheath ratio
r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369         # DMA geometry [m]
Λ₁ = DMAconfig(t, p, qsa, qsa / β, r₁, r₂, l, 0.0, :-, 3, :cylindrical)  # Specify DMA1
Λ₂ = DMAconfig(t, p, qsa, qsa / β, r₁, r₂, l, 0.0, :-, 3, :cylindrical)  # Specify DMA2
bins, z₁, z₂ = 512, dtoz(Λ₁, 500e-9), dtoz(Λ₁, 30e-9) # bins, upper, lower mobility limit
δ₁ = setupDMA(Λ₁, z₁, z₂, bins)                  # Compute matrices
δ₂ = setupDMA(Λ₂, z₁, z₂, bins)                  # Compute matrices

# Upstream Size Distribution
Ax = [[1300.0, 60.0, 1.4], [5000.0, 220.0, 1.6]]
𝕟ᶜⁿ = DMALognormalDistribution(Ax, δ₁)

# Tandem DMA equations
O(k) = mapfoldl(zs -> (δ₂.Ω(Λ₂, δ₂.Z, zs / k, k) .* δ₂.Tl(Λ₂, δ₂.Z, k))', vcat, δ₂.Z)
T₁(zˢ, k) = δ₁.Ω(Λ₁, δ₁.Z, zˢ / k, k) .* δ₁.Tc(k, δ₁.Dp) .* δ₁.Tl(Λ₁, δ₁.Dp)
DMA₁(𝕟, zˢ, gf) = @_ map((gf ⋅ (T₁(zˢ, _) * 𝕟)), 1:3)
itp(𝕟) = interpolateSizeDistributionOntoδ((𝕟, δ₂))
DMA₂(𝕟, k) = O(k) * 𝕟

Dd = 100e-9             # Dry diameter
zˢ = dtoz(Λ₁, Dd);      # Mobility of 100 nm particle
gf = 1.6                # Growth factor
ℕ = DMA₁(𝕟ᶜⁿ, zˢ, gf)   # Transmission through DMA1
𝕄 = map(k -> (@> itp(ℕ[k]) DMA₂(k)), 1:3) # Transmission through DMA2
𝕞ᵗ = sum(𝕄)                               # total response
#hide
mdf(k) = DataFrame(#hide
    Dp = 𝕄[k].Dp./(Dd*1e9), #hide
    S = 𝕄[k].S, #hide
    Dist = ["𝕄[$k]" for i = 1:length(𝕄[k].Dp)]#hide
)#hide
#hide
df1 = mapreduce(mdf, vcat, 1:3)#hide
df2 = DataFrame(Dp = 𝕞ᵗ.Dp./(Dd*1e9), S = 𝕞ᵗ.S, Dist = ["𝕞ᵗ" for i = 1:length(𝕞ᵗ.Dp)])#hide
df = [df2; df1]#hide
#hide
colors = ["black", "darkred", "steelblue3", "darkgoldenrod"]#hide
#hide
p2 = plot(#hide
    df,#hide
    x = :Dp,#hide
    y = :S,#hide
    color = :Dist,#hide
    Geom.line,#hide
    Guide.xlabel("Apparent Growth Factor", orientation = :horizontal),#hide
    Guide.ylabel("dN/dlnD (cm⁻³)"),#hide
    Guide.xticks(ticks = [1.2,1.4,1.6,1.8,2.0]),#hide
    Guide.colorkey(; title = ""),#hide
    Scale.color_discrete_manual(colors...),#hide
    Coord.cartesian(xmin = 1.2, xmax = 2),#hide
    Theme(plot_padding = [5mm, 10mm, 0mm, 0mm]),#hide
)#hide
```

The figure demonstrates the apparent shift toward smaller growth factors for multicharge 
particles. See [Petters (2021)](https://amt.copernicus.org/articles/14/7909/2021/amt-14-7909-2021.pdf) for more explanation. The complete example is reproduced as ```transmission3.jl``` in the ```examples/``` folder of the main repository.

#### Multiple Compositions
The above example can be extended to write a TDMA model that integrates over a pdf. This
function can be obtained from [TDMA1Dpdf](@ref), which is part of the package.

```julia
function TDMA1Dpdf(𝕟ᵢₙ, Λ₁ᵢₙ, Λ₂ᵢₙ, dma2rangeᵢₙ)
    Λ₁, Λ₂, 𝕟1 = deepcopy(Λ₁ᵢₙ), deepcopy(Λ₂ᵢₙ), deepcopy(𝕟ᵢₙ)
    r = deepcopy(dma2rangeᵢₙ)
    Dd, gmin, gmax, n = r[1], r[2], r[3], r[4]
    nDMA, Dmin, Dmax = length(𝕟1.Dp), minimum(𝕟1.Dp), maximum(𝕟1.Dp)

    δ₁ = setupDMA(Λ₁, dtoz(Λ₁, Dmax * 1e-9), dtoz(Λ₁, Dmin * 1e-9), nDMA)
    δ₂ = setupDMA(Λ₂, dtoz(Λ₂, gmax * Dd), dtoz(Λ₂, gmin * Dd), n)
    𝕟 = interpolateSizeDistributionOntoδ((𝕟1, δ₁))

    @memoize O(k) = (hcat(map(i -> δ₂.Ω(Λ₂, δ₂.Z, i/k, k) .* δ₂.Tl(Λ₂, δ₂.Dp), δ₂.Z)...))'
    @memoize T₁(zˢ, k) = δ₁.Ω(Λ₁, δ₁.Z, zˢ / k, k) .* δ₁.Tc(k, δ₁.Dp) .* δ₁.Tl(Λ₁, δ₁.Dp)
	@memoize DMA₁(𝕟, zˢ, gf) = @_ map((gf ⋅ (T₁(zˢ, _) * 𝕟)), 1:6)
	@memoize DMA₂(𝕟, k) = O(k) * 𝕟
	@memoize itp(𝕟) = interpolateSizeDistributionOntoδ((𝕟, δ₂))
	@memoize function TDMA(𝕟, zˢ, gf)
		ℕ = DMA₁(𝕟, zˢ, gf)
		map(k -> (@> itp(ℕ[k]) DMA₂(k)), 1:length(ℕ)) |> sum
	end
	
	@memoize model(𝕟, P, Dd, gf) =
		sum(@_ map(P[_] * TDMA(𝕟, dtoz(Λ₁, Dd), gf[_]), 1:length(P)))
end
```

Note that the basic principle is the same as the single composition above. However,
```DMA₁(𝕟, zˢ, gf)``` sums directly over all charges, so the individual charge distributions are not considered. The function ```TDMA(𝕟, zˢ, gf)``` returns the output from the TDMA. The function ```model(𝕟, P, Dd, gf)``` extends this over a pdf, where gf is a list of growth fractors and P are corresponding probabilities. It one possible implementation of Eqs. (16) and (17) in [Petters (2021)](https://amt.copernicus.org/articles/14/7909/2021/amt-14-7909-2021.pdf)

Below is an example with 4 population each having a unique growth factor and fractional contribution to the total distribution. If the fractions are known, the net response function of the TDMA is readily computed. The example is reproduced as ```transmission4.jl``` in the examples folder of the main repository.  

```@example
using Distributions #hide
using DifferentialMobilityAnalyzers #hide
using Gadfly #hide
using Printf #hide
t, p = 295.15, 1e5
qsa, qsh = 1.66e-5, 8.33e-5
r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369
Λ₁ = DMAconfig(t, p, qsa, qsh, r₁, r₂, l, 0.0, :-, 6, :cylindrical)
Λ₂ = DMAconfig(t, p, qsa, qsh, r₁, r₂, l, 0.0, :-, 6, :cylindrical)
bins, z₁, z₂ = 120, dtoz(Λ₁, 500e-9), dtoz(Λ₁, 30e-9) # bins, upper, lower mobility limit
δ₁ = setupDMA(Λ₁, z₁, z₂, bins)                

Ax = [[1300.0, 60.0, 1.4], [5000.0, 220.0, 1.6]] 
𝕟 = DMALognormalDistribution(Ax, δ₁)

# scan 100 nm Dd from 0.8Dd to 3.0Dd with 100 bins
dma2range = (100e-9, 0.8, 3.0, 120)

# Get the model function
model = TDMA1Dpdf(𝕟, Λ₁, Λ₂, dma2range)

P = [0.5,0.15, 0.10, 0.25]   # Probability of growth factor (4 populations)
gf = [1.0, 1.2, 1.6, 2.1]    # Values of growth factor
𝕘 = model(𝕟, P, dma2range[1], gf)
#hide
set_default_plot_size(14cm, 8cm)#hide
xlabels = collect(1:0.5:3)#hide
p1 = plot(#hide
    x = 𝕘.Dp./100.0,#hide
    y = 𝕘.N,#hide
    Geom.step,#hide
    Guide.xlabel("Growth Factor (-)"),#hide
    Guide.ylabel("Number concentration (cm-3)", orientation = :vertical),#hide
    Guide.xticks(ticks = (collect(0.8:0.1:3))),#hide
    Scale.x_continuous(labels = x -> x in xlabels ? @sprintf("%.1f", (x)) : ""),#hide
    Coord.cartesian(xmin = 0.8, xmax = 3),#hide
    Theme(plot_padding = [2mm, 2mm, 2mm, 2mm]),#hide
)#hide
```


### Volatilty Tandem DMA

See the [Manuscript](https://www.tandfonline.com/doi/full/10.1080/02786826.2018.1530724) and notebook S9 in the [Notebooks](@ref) section for examples involving the volatility tandem DMA. 

### Dual Tandem DMA

See the [Manuscript](https://www.tandfonline.com/doi/full/10.1080/02786826.2018.1530724) and notebook S10 and S11 in the [Notebooks](@ref) section for examples involving the volatility tandem DMA. 

