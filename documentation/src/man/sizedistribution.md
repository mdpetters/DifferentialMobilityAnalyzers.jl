# Size Distributions

## Notation
Since this package is working a lot with size distributions, it is useful to highlight them. By convention we will use blackboard bold characters to denote size distributions

!!! info
    **Blackboard bold font**: 𝕒, 𝕓, 𝕔, 𝕕, 𝕖, 𝕗, 𝕘, ...

## Histogram Representation
The size distribution is represented as a histogram. It is a composite data type [SizeDistribution](@ref) that has the fields Dp, De, ΔlnD, S, and N. 

- 𝕟.Dp: Geometric midpoint diameters
- 𝕟.De: Bin edge diameters
- 𝕟.ΔlnD: Log bin spacing, ΔlnD = ln(Dup/Dlow)
- 𝕟.S: Spectral density
- 𝕟.N: Number concentration in the bin

The size distribution can be constructed various ways. The easiest is to use one of the constructor functions. For example, the [lognormal](@ref) function creates a lognormal size distribution. The following example creates a lognormal size distribution with number concentration equals 200 cm-3, geometric mode diameter of 80 nm, geometric standard deviation of 1.2, with 10 size bins between 30 and 300 nm. The result is placed in a DataFrame for display purposes. The ```r``` function is to round the results for clarity. The output illustrates that the SizeDistribution type is simply a histagram table.

```@example
using DifferentialMobilityAnalyzers, DataFrames # hide
r(x) = round.(Int,x)    # Function to round and convert to Int
𝕟 = lognormal([[200, 80, 1.2]]; d1 = 30.0, d2 = 300.0, bins = 10);
DataFrame(
    Dlow = r(𝕟.De[1:end-1]),
    Dup = r(𝕟.De[2:end]),
    ΔlnD = round.(𝕟.ΔlnD, digits = 2),
    Dp = r(𝕟.Dp),
    S = r(𝕟.S),
    N = r(𝕟.N),
)
```

## Manipulating Size Distributions
Size distributions can be intuitively manipulated through operators. For example, the sum of two size distributions (𝕩 = 𝕟₁ + 𝕟₂) is the superposition. 

```@example
using DifferentialMobilityAnalyzers, Gadfly, DataFrames, Printf # hide
# Example addition of size distributions
𝕟₁ = lognormal([[120, 90, 1.20]]; d1 = 10.0, d2 = 1000.0, bins = 256)   # size distribution
𝕟₂ = lognormal([[90, 140, 1.15]]; d1 = 20.0, d2 = 800.0, bins = 256)    # size distribution
𝕩 = 𝕟₁ + 𝕟₂
# hide
function getplots1(𝕟₁, 𝕟₂, 𝕩)# hide
    df1 = DataFrame(Dp = 𝕟₁.Dp, S = 𝕟₁.S, Dist = ["𝕟₁" for i = 1:length(𝕟₁.Dp)])# hide
    df2 = DataFrame(Dp = 𝕟₂.Dp, S = 𝕟₂.S, Dist = ["𝕟₂" for i = 1:length(𝕟₂.Dp)])# hide
    df3 = DataFrame(Dp = 𝕩.Dp, S = 𝕩.S, Dist = ["𝕩" for i = 1:length(𝕩.Dp)]) # hide
    df = [df1; df2; df3] # hide
 # hide
    xlabels = log10.([40, 100, 400]) # hide
    colors = ["darkred", "steelblue3", "black"] # hide
# hide
    set_default_plot_size(16cm, 9cm) # hide
    return plot( # hide
        df, # hide
        x = :Dp, # hide
        y = :S, # hide
        color = :Dist, # hide
        Geom.step, # hide
        Guide.xlabel("Particle diameter (nm)"), # hide
        Guide.ylabel("dN/dlnD (cm-3)"), # hide
        Guide.xticks(ticks = log10.([40, 50, 60, 70, 80, 90, 100, 200, 300, 400])), # hide
        Guide.colorkey(; title = "Distribution"), # hide
        Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""), # hide
        Scale.color_discrete_manual(colors...), # hide
        Coord.cartesian(xmin = log10(40), xmax = log10(400)), # hide
    ) # hide
end # hide

getplots1(𝕟₁, 𝕟₂, 𝕩) # hide
```

The package implements a list of [Operators](@ref) for size distribution manipulation. Check out the [Tutorial](@ref) Session 1 and/or Notebook S3 in the [Notebooks](@ref) section for visualizations.