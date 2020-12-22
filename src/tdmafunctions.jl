@doc raw"""
    gfₖ(Λ, zˢ, gf, k)  

Returns the effective growth factor for multi-charge particles. 
- Λ  - DMA configuration 
- zˢ - the dry particle diameter mobility [m2 s-1 V-1]
- gf - the true growth factor
- k  - the number of charges    
    
The effect is described by Gysel et al. (2009) and Shen et al. (2020). 
Shen et al. (2020) write: "For electrical mobility diameter of 100 nm, the doubly and 
triply charged particles are about 151 nm and 196 nm, respectively. When all these three 
kind of particles have a true growth factor of 1.6, they will grow to the size of 160 nm, 
242 nm and 314 nm. Since the number of charges they carry remain the same as before, 
their peak sizes in the second DMA are around 160 nm, 154 nm and 150 nm. Therefore, the 
growth factors they display in the HTDMA measurement is 1.6, 1.54 and 1.5, respectively."

Example use
```
# Define a DMA
t, p = 295.15, 1e5
qsa, qsh = 1.66e-5, 8.33e-5
r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369
Λ = DMAconfig(t, p, qsa, qsh, r₁, r₂, l, 0.0, :-, 6, :cylindrical)

zˢ = dtoz(Λ, 100e-9) # mobility z-star for 100 nm particle 
truegf = 1.6          # true growth factor of 1.6
gfₖ(Λ, zˢ, truegf, 1)    # effective growth factor for 1 charge = 1.6
gfₖ(Λ, zˢ, truegf, 2)    # effective growth factor for 2 charges = 1.544 
gfₖ(Λ, zˢ, truegf, 3)    # effective growth factor for 3 charges = 1.507 
gfₖ(Λ, zˢ, truegf, 4)    # effective growth factor for 4 charges = 1.481 
```
"""
@memoize gfₖ(Λ, zˢ, gf, k) = ztod(Λ, 1, dtoz(Λ, 1e-9 * ztod(Λ, k, zˢ) * gf) * k) ./ ztod(Λ, 1, zˢ) 

@doc raw"""
    TDMA1Dpdf(𝕟ᵢₙ,  Λ₁ᵢₙ , Λ₂ᵢₙ, dma2rangeᵢₙ)

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
bins, z₁, z₂ = 120, dtoz(Λ₁, 500e-9), dtoz(Λ₁, 30e-9) # bins, upper, lower mobility limit
δ₁ = setupDMA(Λ₁, z₁, z₂, bins)                
Ax = [[1300.0, 60.0, 1.4], [5000.0, 220.0, 1.6]] 
𝕟 = DMALognormalDistribution(Ax, δ₁)

# scan 100 nm Dd from 0.8Dd to 3.0Dd with 100 bins
dma2range = (100e-9, 0.8, 3.0, 30)

# Get the model function
model = TDMA1Dpdf(𝕟, Λ₁, Λ₂, dma2range)

# Growth factor grid along with the PDF is evaluated over
Ax = [[1300.0, 60.0, 1.4], [5000.0, 220.0, 1.6]] 
𝕟 = DMALognormalDistribution(Ax, δ₁)

P = [0.5,0.15, 0.10, 0.25]   # Probability of growth factor (4 populations)
gf = [1.0, 1.2, 1.6, 2.1]    # Values of growth factor
𝕘 = mymodel(𝕟, P, dma2range[1], gf)

plot(x = 𝕘.Dp/(dma2range[1]*1e9), y = 𝕘.N, Geom.line,
    Guide.xticks(ticks = 0.8:0.2:3),
    Coord.cartesian(xmin = 0.8, xmax = 3.0))
```
"""
function TDMA1Dpdf(𝕟ᵢₙ,  Λ₁ᵢₙ , Λ₂ᵢₙ, dma2rangeᵢₙ)
    Λ₁ , Λ₂, 𝕟1 = deepcopy(Λ₁ᵢₙ), deepcopy(Λ₂ᵢₙ), deepcopy(𝕟ᵢₙ)
    r = deepcopy(dma2rangeᵢₙ)
    Dd, gmin, gmax, n = r[1], r[2], r[3], r[4]
    nDMA, Dmin, Dmax = length(𝕟1.Dp), minimum(𝕟1.Dp), maximum(𝕟1.Dp)

    δ₁ = setupDMA(Λ₁, dtoz(Λ₁, Dmax*1e-9), dtoz(Λ₁, Dmin*1e-9), nDMA)
	δ₂ = setupDMA(Λ₂, dtoz(Λ₂, gmax*Dd), dtoz(Λ₂, gmin*Dd), n)
    𝕟 = interpolateSizeDistributionOntoδ((𝕟1, δ₁))
    
	@memoize T₁(zˢ, k) = δ₁.Ω(Λ₁, δ₁.Z, zˢ / k) .* δ₁.Tc(k, δ₁.Dp) .* δ₁.Tl(Λ₁, δ₁.Dp)
    @memoize cr(zˢ, k) = ztod(Λ₁, 1, zˢ) / ztod(Λ₁, k, zˢ)
    @memoize DMA₁(𝕟, zˢ, gf) = sum(@_ map(cr(zˢ, _) ⋅ (gfₖ(Λ₁, zˢ, gf, _) ⋅ (T₁(zˢ, _) * 𝕟)), 1:6))
    @memoize DMA₂(𝕟) = δ₂.𝐎 * 𝕟
	@memoize itp(𝕟) = interpolateSizeDistributionOntoδ((𝕟, δ₂))
	@memoize TDMA(𝕟, zˢ, gf) = @> DMA₁(𝕟, zˢ, gf) itp DMA₂
	@memoize model(𝕟, P, Dd, gf) = sum(@_ map(P[_]*TDMA(𝕟, dtoz(Λ₁, Dd), gf[_]), 1:length(P)))
end

@doc raw"""
    TDMA1Ddomainfunction(𝕟ᵢₙ,  Λ₁ᵢₙ , Λ₂ᵢₙ, dma2rangeᵢₙ)

Returns a domain function that can be used with RegularizationTools.jl to create
a [domain matrix](https://mdpetters.github.io/RegularizationTools.jl/stable/manual/#Creating-a-Design-Matrix).
The inputs are same as for [TDMA1Dpdf](@ref) to compute a model that models the output of a tandem DMA for an 
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
model = TDMA1Dpdf(𝕟, Λ₁, Λ₂, dma2range)

# Growth factor grid along with the PDF is evaluated over
mgf = 0.8:0.05:2.5 

# A growth factor PDF
gfpdf = pdf(truncated(Normal(1.2,0.2) , 1, 17), mgf)

# The output model
𝕘 = model(𝕟, gfpdf, dma2range[1], mgf)

# The domain function
f = TDMA1Ddomainfunction(𝕟, Λ₁, Λ₂, dma2range)

# Computes the design matrix
A = designmatrix(mgf, f)

# Computes the output via matrix multiplication
g = A*gfpdf

# Compare with the regular model. Note that the matrix method follows the mgf grid.
# which can be different than the model grid given via dma2range
plot(layer(x = 𝕘.Dp/(dma2range[1]*1e9), y = 𝕘.N), layer(x = mgf, y = g, Geom.line))
```
"""
function TDMA1Ddomainfunction(𝕟ᵢₙ,  Λ₁ᵢₙ , Λ₂ᵢₙ, dma2rangeᵢₙ)
    Λ₁ , Λ₂ = deepcopy(Λ₁ᵢₙ), deepcopy(Λ₂ᵢₙ)
    𝕟 = deepcopy(𝕟ᵢₙ)
    r = deepcopy(dma2rangeᵢₙ)
    Dd = r[1]
    model = TDMA1Dpdf(𝕟ᵢₙ,  Λ₁ᵢₙ , Λ₂ᵢₙ, dma2rangeᵢₙ)

    function f(domain::Domain)
        gf, P, ogf = domain.s, domain.x, domain.q
        out = model(𝕟, P, Dd, gf)
        mgf = reverse(out.Dp)./(Dd*1e9)
        N = reverse(out.N)
        fitp = @> interpolate((mgf,), N, Gridded(Linear())) extrapolate(0.0)
        return fitp(ogf)  
    end
end
