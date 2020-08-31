# Transmission Through the DMA at a Constant Voltage

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
    Guide.xlabel("Particle diameter (nm)"), # hide
    Guide.ylabel("dN/dlnD (cm-3)"), # hide
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