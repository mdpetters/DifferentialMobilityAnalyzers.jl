@doc raw"""
    TDMA_1D_pdf(𝕟ᵢₙ,  Λ₁ᵢₙ , Λ₂ᵢₙ, dma2rangeᵢₙ)

Returns a function model that models the output of a tandem DMA for an 
input pdf of growth factors. Tne function is specialized for a specific size distribution

- 𝕟ᵢₙ is the size distribution with  Dp sorted in in ascending order and units of nm
- Λ₁ᵢₙ , Λ₂ᵢₙ DMA 1 and 2 configuration of type of [DMAconfig](@ref)
- dma2rangeᵢₙ a tuple (Dd, gmin, gmax, n) where
    Dd is the dry diameter selected by DMA 1 in units of m
    gmin is the lower range of growth factors scanned by DMA2
    gmax is the upper range of growth factor scanned by DMA2
    n is the number of bins to represent the DMA2 grid

Example use
```julia
using Distributions
using DifferentialMobilityAnalyzers
using Gadfly

t, p = 295.15, 1e5
qsa, qsh = 1.66e-5, 8.33e-5
r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369
Λ₁ = DMAconfig(t, p, qsa, qsh, r₁, r₂, l, 0.0, :-, 6, :cylindrical)
Λ₂ = DMAconfig(t, p, qsa, qsh, r₁, r₂, l, 0.0, :-, 6, :cylindrical)
𝕟 = lognormal([[9., 40., 1.5], [500., 180., 1.4]]; d1 = 800, d2 = 10.0, bins = 120)

# scan 100 nm Dd from 0.8Dd to 3.0Dd with 100 bins
dma2range = (100e-9, 0.8, 3.0, 30)

# Get the model function
model = TDMA_1D_pdf(𝕟, Λ₁, Λ₂, dma2range)

# Growth factor grid along with the PDF is evaluated over
mgf = 0.8:0.05:2.5 

# A growth factor PDF
gfpdf = pdf(truncated(Normal(1.2,0.2) , 1, 17), mgf)

# The output model
𝕘 = model(𝕟, gfpdf, dma2range[1], mgf)


plot(x = 𝕘.Dp/(dma2range[1]*1e9), y = 𝕘.N)
```
"""
function TDMA_1D_pdf(𝕟ᵢₙ,  Λ₁ᵢₙ , Λ₂ᵢₙ, dma2rangeᵢₙ)
    Λ₁ , Λ₂, 𝕟1 = deepcopy(Λ₁ᵢₙ), deepcopy(Λ₂ᵢₙ), deepcopy(𝕟ᵢₙ)
    r = deepcopy(dma2rangeᵢₙ)
    Dd, gmin, gmax, n = r[1], r[2], r[3], r[4]
    nDMA, Dmin, Dmax = length(𝕟1.Dp), minimum(𝕟1.Dp), maximum(𝕟1.Dp)

    δ₁ = setupDMA(Λ₁, dtoz(Λ₁, Dmax*1e-9), dtoz(Λ₁, Dmin*1e-9), nDMA)
	δ₂ = setupDMA(Λ₂, dtoz(Λ₂, gmax*Dd), dtoz(Λ₂, gmin*Dd), n)
    𝕟 = interpolateSizeDistributionOntoδ((𝕟1, δ₁))
    
	@memoize T₁(zˢ, k) = δ₁.Ω(Λ₁, δ₁.Z, zˢ / k) .* δ₁.Tc(k, δ₁.Dp) .* δ₁.Tl(Λ₁, δ₁.Dp)
	@memoize cr(k, zˢ) = ztod(Λ₁, 1, zˢ) / ztod(Λ₁, k, zˢ)
	@memoize DMA₁(𝕟, zˢ) = sum(@_ map(cr(_, zˢ) ⋅ (T₁(zˢ, _) * 𝕟), 1:6))
	@memoize DMA₂(𝕟, gf) = δ₂.𝐎 * (gf ⋅ 𝕟)
	@memoize itp(𝕟) = interpolateSizeDistributionOntoδ((𝕟, δ₂))
	@memoize TDMA(𝕟, zˢ, gf) = @> DMA₁(𝕟, zˢ) itp DMA₂(gf)
	@memoize model(𝕟, P, Dd, gf) = sum(@_ map(P[_]*TDMA(𝕟, dtoz(Λ₁, Dd), gf[_]), 1:length(P)))
end

@doc raw"""
    TDMA_1D_domainfunction(𝕟ᵢₙ,  Λ₁ᵢₙ , Λ₂ᵢₙ, dma2rangeᵢₙ)

Returns a domain function that can be used with RegularizationTools.jl to create
a [domain matrix](https://mdpetters.github.io/RegularizationTools.jl/stable/manual/#Creating-a-Design-Matrix).
The inputs are same as for [TDMA_1D_pdf](@ref) to compute a model that models the output of a tandem DMA for an 
input pdf of growth factors. Tne function and matrix are specialized for a specific size distribution
    
- 𝕟ᵢₙ is the size distribution with  Dp sorted in in ascending order and units of nm
- Λ₁ᵢₙ , Λ₂ᵢₙ DMA 1 and 2 configuration of type of [DMAconfig](@ref)
- dma2rangeᵢₙ a tuple (Dd, gmin, gmax, n) where
    Dd is the dry diameter selected by DMA 1 in units of m
    gmin is the lower range of growth factors scanned by DMA2
    gmax is the upper range of growth factor scanned by DMA2
    n is the number of bins to represent the DMA2 grid

Example
```julia
using Distributions
using DifferentialMobilityAnalyzers
using RegularizationTools
using Gadfly

t, p = 295.15, 1e5
qsa, qsh = 1.66e-5, 8.33e-5
r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369
Λ₁ = DMAconfig(t, p, qsa, qsh, r₁, r₂, l, 0.0, :-, 6, :cylindrical)
Λ₂ = DMAconfig(t, p, qsa, qsh, r₁, r₂, l, 0.0, :-, 6, :cylindrical)
𝕟 = lognormal([[9., 40., 1.5], [500., 180., 1.4]]; d1 = 800, d2 = 10.0, bins = 120)

# scan 100 nm Dd from 0.8Dd to 3.0Dd with 100 bins
dma2range = (100e-9, 0.8, 3.0, 30)

# Get the model function
model = TDMA_1D_pdf(𝕟, Λ₁, Λ₂, dma2range)

# Growth factor grid along with the PDF is evaluated over
mgf = 0.8:0.05:2.5 

# A growth factor PDF
gfpdf = pdf(truncated(Normal(1.2,0.2) , 1, 17), mgf)

# The output model
𝕘 = model(𝕟, gfpdf, dma2range[1], mgf)

# The domain function
f = TDMA_1D_domainfunction(𝕟, Λ₁, Λ₂, dma2range)

# Computes the design matrix
A = designmatrix(mgf, f)

# Computes the output via matrix multiplication
g = A*gfpdf

# Compare with the regular model. Note that the matrix method follows the mgf grid.
# which can be different than the model grid given via dma2range
plot(layer(x = 𝕘.Dp/(dma2range[1]*1e9), y = 𝕘.N), layer(x = mgf, y = g, Geom.line))
```
"""
function TDMA_1D_domainfunction(𝕟ᵢₙ,  Λ₁ᵢₙ , Λ₂ᵢₙ, dma2rangeᵢₙ)
    Λ₁ , Λ₂ = deepcopy(Λ₁ᵢₙ), deepcopy(Λ₂ᵢₙ)
    𝕟 = deepcopy(𝕟ᵢₙ)
    r = deepcopy(dma2rangeᵢₙ)
    Dd = r[1]
    model = TDMA_1D_pdf(𝕟ᵢₙ,  Λ₁ᵢₙ , Λ₂ᵢₙ, dma2rangeᵢₙ)

    function f(domain::Domain)
        gf, P, ogf = domain.s, domain.x, domain.q
        out = model(𝕟, P, Dd, gf)
        mgf = reverse(out.Dp)./(Dd*1e9)
        N = reverse(out.N)
        fitp = @> interpolate((mgf,), N, Gridded(Linear())) extrapolate(0.0)
        return fitp(ogf)  
    end
end
