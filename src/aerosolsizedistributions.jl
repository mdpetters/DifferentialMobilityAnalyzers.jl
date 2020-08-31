# +
# This file handles binned size distribution from size distribution functions
#
# Author: Markus Petters (mdpetter@ncsu.edu)
# 	      Department of Marine Earth and Atmospheric Sciences
#         NC State University
#         Raleigh, NC 27605
#
#         April, 2018
#-


# Helper functions for lognormal distribution
md(A, x) = @. A[1] / (√(2π) * log(A[3])) * exp(-(log(x / A[2]))^2 / (2log(A[3])^2))
logn(A, x) = mapreduce((A) -> md(A, x), +, A)

@doc raw"""
    lognormal(A; d1 = 8.0, d2 = 2000.0, bins = 256)
    

The lognormal function instantiates a the [SizeDistribution](@ref) type with a multi-modal
lognormal distribution. The multi-modal lognormal size distribution is given by 
(e.g. Seinfeld and Pandis, 2006)

``\frac{dN}{d\ln D_p} = \sum_{i=1}^n \frac{N_{t,i}}{\sqrt{2\pi}\ln\sigma_{g,i}} 
\exp \left(- \frac{\left[\ln D_p-\ln D_{pg,i}\right]^2}{2\ln \sigma_{g,i}^2}\right)``

where ``\frac{dN}{d\ln D_p}`` is the spectral number density, ``N_{t,i}`` is the total
number concentration, ``\sigma_{g,i}`` is the geometric standard deviation,  
``D_{pg,i}`` is the geometric mean diameter of the ``i^{th}`` mode,  ``n``is the number of 
modes.  
    
Each mode is coded as an array of [Nt, Dg, sg]. The inputs are
- A is an array of arrays with modes, i.e. [[Nt1,Dg1,sg1],[Nt2,Dg2,sg2], ...]
- d1 is the lower diameter of the grid
- d2 is the upper diameter of the grid
- bins is the number of size bins.

By definition of the function sg >= 1, with sg1 corresponding to an infinitely narrow mode
The function is unit agnostic. 

Example Usage
```julia
𝕟 = lognormal([[200.0, 80.0, 1.3]]; d1 = 10, d2 = 500.0, bins = 120)
𝕟 = lognormal([[200.0, 80.0, 1.3], [200.0, 150.0, 1.3]]; d1 = 10, d2 = 800.0, bins = 60)
```
"""
function lognormal(A; d1 = 8.0, d2 = 2000.0, bins = 256)
    De = 10.0 .^ range(log10(d1), stop = log10(d2), length = bins + 1)
    Dp = sqrt.(De[2:end] .* De[1:end-1])
    ΔlnD = log.(De[2:end] ./ De[1:end-1])
    S = logn(A, Dp)
    N = S .* ΔlnD
    return SizeDistribution(A, De, Dp, ΔlnD, S, N, :lognormal)
end

@doc raw"""
    triangular(Λ::DMAconfig, δ::DifferentialMobilityAnalyzer, A)

Instantiates a single mode triangular distribution in mobility space with number 
concentration Nt and mode diameter Dg. This is a convenient constructor to model a single
mode of the distribution output of an idealized DMA.

Example Usage
```julia
𝕟 = triangular(Λ, δ, [200.0, 50.0])
```
"""
function triangular(Λ::DMAconfig, δ::DifferentialMobilityAnalyzer, A)
    Nt = A[1]
    zˢ = dtoz(Λ, A[2] .* 1e-9)
    Ntrans = Ω(Λ, δ.Z, zˢ)  # Transfer function with peak at 1
    N = Nt .* Ntrans ./ sum(Ntrans)
    S = N ./ δ.ΔlnD
    return SizeDistribution(
        [A],
        reverse(δ.De),
        reverse(δ.Dp),
        reverse(δ.ΔlnD),
        reverse(S),
        reverse(N),
        :DMA,
    )
end

@doc raw"""
    DMALognormalDistribution(A, δ::DifferentialMobilityAnalyzer)
    
The DMALognormalDistribution function instantiates a the [SizeDistribution](@ref) type with
a multi-modal lognormal distribution. The multi-modal lognormal size distribution is given by 
(e.g. Seinfeld and Pandis, 2006)

``\frac{dN}{d\ln D_p} = \sum_{i=1}^n \frac{N_{t,i}}{\sqrt{2\pi}\ln\sigma_{g,i}} 
\exp \left(- \frac{\left[\ln D_p-\ln D_{pg,i}\right]^2}{2\ln \sigma_{g,i}^2}\right)``

where ``\frac{dN}{d\ln D_p}`` is the spectral number density, ``N_{t,i}`` is the total
number concentration, ``\sigma_{g,i}`` is the geometric standard deviation,  
``D_{pg,i}`` is the geometric mean diameter of the ``i^{th}`` mode,  ``n``is the number of 
modes.  
    
Each mode is coded as an array of [Nt, Dg, sg]. The inputs are
- A is an array of arrays with modes, i.e. [[Nt1,Dg1,sg1],[Nt2,Dg2,sg2], ...]
- δ is a DifferentialMobilityAnalyzer

By definition of the function sg >= 1, with sg1 corresponding to an infinitely narrow mode
The function is unit agnostic. The diameter grid is that from the 
DifferentialMobilityAnalyzer δ. 

Example Usage
```julia
t,p = 295.15, 1e5                        
qsa,qsh = 1.66e-5, 8.3e-5                     
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369               
Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,0.0,:-,6,:cylindrical) 
bins,z₁,z₂ = 30, vtoz(Λ,10000), vtoz(Λ,10)   
δ = setupDMA(Λ, z₁, z₂, bins)

𝕟 = DMALognormalDistribution([[200.0, 80.0, 1.3]], δ)
𝕟 = DMALognormalDistribution([[200.0, 80.0, 1.3], [200.0, 150.0, 1.3]], δ)
```
"""
function DMALognormalDistribution(A, δ::DifferentialMobilityAnalyzer)
    S = logn(A, δ.Dp)

    return SizeDistribution(A, δ.De, δ.Dp, δ.ΔlnD, S, S .* δ.ΔlnD, :DMA)
end

"""
    *(a::Number, 𝕟::SizeDistribution)

Multiplication of scalar and size distribution. The net result is a scaling of the 
number concentration of the spectra by a. The function is symmetric such that a * 𝕟 == 𝕟 * a.

Let a denote a number and 𝕟 denote a size distribution. Then
```julia
𝕩 = a * 𝕟 
```
is defined such that

```julia
𝕩.N = a * 𝕟.N
𝕩.S = a * 𝕟.S
```

Example Usage
```julia
𝕟 = lognormal([[120, 90, 1.20]]; d1 = 10.0, d2 = 1000.0, bins = 256)
𝕩 = 2.3 * 𝕟₂
```
"""
function *(a::Number, 𝕟::SizeDistribution)
    # This function defines the product of a scalar and a size distribution
    N = a * 𝕟.N
    S = a * 𝕟.S

    return SizeDistribution([[]], 𝕟.De, 𝕟.Dp, 𝕟.ΔlnD, S, N, :axdist)
end

*(𝕟::SizeDistribution, a::AbstractFloat) = *(a::AbstractFloat, 𝕟::SizeDistribution)


"""
    *(a::Vector{<:AbstractFloat}, 𝕟::SizeDistribution)

Multiplication of vector and size distribution. The net result is a bin-by-bin scaling 
of the number concentration. The function is symmetric such that a * 𝕟 == 𝕟 * a.

Let T denote a 1D vector that has the same number of elements as the 
size distribution 𝕟. Then
```julia
𝕩 = T * 𝕟 
```
is defined such that 
```julia
𝕩.N = T * 𝕟.N
𝕩.S = T * 𝕟.S
```
    
Example Usage
```julia
𝕟 = lognormal([[100, 100, 1.1]]; d1 = 10.0, d2 = 1000.0, bins = 256)  
μ,σ = 100.0, 200.0
T = 0.5*(1.0 .+ erf.((𝕟.Dp .- μ)./(sqrt(2σ)))) 
𝕩 = T * 𝕟                                        
```
"""
function *(a::Vector{<:AbstractFloat}, 𝕟::SizeDistribution)
    # This function defines the product of a vector and a size distribution
    N = a .* 𝕟.N
    S = a .* 𝕟.S
    return SizeDistribution([[]], 𝕟.De, 𝕟.Dp, 𝕟.ΔlnD, S, N, :axdist)
end

*(𝕟::SizeDistribution, a::Vector{<:AbstractFloat}) =
    *(a::Vector{<:AbstractFloat}, 𝕟::SizeDistribution)


"""
    *(𝐀::AbstractMatrix, 𝕟::SizeDistribution)

Multiplication of matrix and size distribution. The net result is the multiplication of the 
matrix with number concentration and spectral density fields. 

Let 𝐀 denote an nxn matrix where n equals the number of size bins of 𝕟. Then
```julia
𝕩 = 𝐀 * 𝕟 
```
is defined such that 
```julia
𝕩.N = 𝐀 * 𝕟.N
𝕩.S = 𝐀 * 𝕟.S
```

```julia
𝕟 = lognormal([[100, 100, 1.1]]; d1 = 10.0, d2 = 1000.0, bins = 256)
𝐀 = rand(256,256)
𝕩 = 𝐀 * 𝕟                                        
```
"""
function *(𝐀::AbstractMatrix, 𝕟::SizeDistribution)
    # This function defines the product of a matrix and a size distribution
    N = 𝐀 * 𝕟.N
    S = 𝐀 * 𝕟.S
    return SizeDistribution([[]], 𝕟.De, 𝕟.Dp, 𝕟.ΔlnD, S, N, :axdist)
end


"""
    *(𝕟₁::SizeDistribution, 𝕟₂::SizeDistribution)

Multiplication of size distribution and a size distribution. The net result is a 
size distribution that has total number concentration square. For a probability 
distributions that by definition integrate to unity, this operation corresponds to 
the product of two random variates with distribution 1 and 2.
    
Let 𝕟₁ and 𝕟₂ denote a two size distribution defined on the same diameter grid. Then
```julia
𝕩 = 𝕟₁ * 𝕟₂ 
```
is defined such that 
```julia
Nsq = 𝕟₁.N * 𝕟₂.N
𝕩.N = sum(𝕟₁.N) * sum(𝕟₂.N) * Nsq./sum(Nsq)
𝕩.S = N ./ 𝕟₁.ΔdlnD
```

Example Usage
```julia
𝕟₁ = lognormal([[120, 90, 1.20]]; d1 = 10.0, d2 = 1000.0, bins = 256)
𝕟₂ = lognormal([[90, 140, 1.15]]; d1 = 20.0, d2 = 800.0, bins = 64)
𝕩 = 𝕟₁ * 𝕟₂
```
"""
function *(𝕟₁::SizeDistribution, 𝕟₂::SizeDistribution)
    # This function defines the product of two size distributions
    Nsq = 𝕟₁.N .* 𝕟₂.N
    N = sum(𝕟₁.N) * sum(𝕟₂.N) * Nsq ./ sum(Nsq)
    S = N ./ 𝕟₁.ΔlnD
    return SizeDistribution([[]], 𝕟₁.De, 𝕟₁.Dp, 𝕟₁.ΔlnD, S, N, :dist_sq)
end

"""
    /(𝕟₁::SizeDistribution, 𝕟₂::SizeDistribution)

Division of size distribution and size distribution. The net result is a size distribution 
that is the ratio of the concentration vectors.

Let 𝕟₁ and 𝕟₂ denote a two size distribution defined on the same diameter grid. Then
    
```julia
𝕩 = 𝕟₁ / 𝕟₂
```
is defined such that

```julia
N = 𝕟₁.N ./ 𝕟₂.N
S = 𝕟₁.S ./ 𝕟₂.S
```

Example Usage
```julia
𝕟₁ = lognormal([[120, 90, 1.20]]; d1 = 10.0, d2 = 1000.0, bins = 256)
𝕟₂ = lognormal([[90, 140, 1.15]]; d1 = 20.0, d2 = 800.0, bins = 64)
𝕩 = 𝕟₁ / 𝕟₂
```

"""
function /(𝕟₁::SizeDistribution, 𝕟₂::SizeDistribution)
    # This function defines the product of two size distributions
    N = 𝕟₁.N ./ 𝕟₂.N
    S = 𝕟₁.S ./ 𝕟₂.S
    return SizeDistribution([[]], 𝕟₁.De, 𝕟₁.Dp, 𝕟₁.ΔlnD, S, N, :dist_sq)
end


"""
    ⋅(a::Number, 𝕟::SizeDistribution)

Multiplication of a scalar and a size distribution. The net result is a uniform diameter 
shift of the size distribution. The function is symmetric such that a ⋅ 𝕟 == 𝕟 ⋅ a.

Let a denote a floating point scalar and 𝕟 denote a size distribution. Then
```julia
𝕩 = a ⋅ 𝕟
```
is defined such that 
```julia
𝕩.Dp = a * 𝕟.Dp 
```

Example Usage
```julia
a = 2.0 
𝕟 = lognormal([[300, 100, 1.3]]; d1 = 10.0, d2 = 1000.0, bins = 256)
𝕩 = a ⋅ 𝕟 
```
"""
function LinearAlgebra.:⋅(a::Number, 𝕟::SizeDistribution)
    if 𝕟.Dp[1] > 𝕟.Dp[2]
        nDp = reverse(a * 𝕟.Dp)
        itpN = interpolate((nDp,), reverse(𝕟.N), Gridded(Linear()))
        extN = extrapolate(itpN, 0)

        itpS = interpolate((nDp,), reverse(𝕟.S), Gridded(Linear()))
        extS = extrapolate(itpS, 0)
        N = clean(extN(reverse(𝕟.Dp)))
        S = clean(extS(reverse(𝕟.Dp)))
        N = S .* reverse(𝕟.ΔlnD)
        return SizeDistribution([[]], 𝕟.De, 𝕟.Dp, 𝕟.ΔlnD, reverse(S), reverse(N), :axdist)
    else
        nDp = a * 𝕟.Dp
        itpN = interpolate((nDp,), reverse(𝕟.N), Gridded(Linear()))
        extN = extrapolate(itpN, 0)

        itpS = interpolate((nDp,), reverse(𝕟.S), Gridded(Linear()))
        extS = extrapolate(itpS, 0)
        N = clean(extN(𝕟.Dp))
        S = clean(extS(𝕟.Dp))

        N = S .* 𝕟.ΔlnD
        return SizeDistribution([[]], 𝕟.De, 𝕟.Dp, 𝕟.ΔlnD, S, N, :axdist)
    end
end

⋅(𝕟::SizeDistribution, a::AbstractFloat) = ⋅(a::AbstractFloat, 𝕟::SizeDistribution)

"""
    LinearAlgebra.:⋅(A::Vector{<:AbstractFloat}, 𝕟::SizeDistribution)

Dot product of vector and size distribution.   The net result is diameter dependent shift of 
the size distribution. The function is symmetric such that A ⋅ 𝕟 == 𝕟 ⋅ A.


Let T denote a floating point vector with the same number of elements as the size distribution 𝕟. Then
```julia
𝕩 = T ⋅ 𝕟 
```

is defined such that 
```julia
𝕩.Dp = T .* 𝕟.dp 
```

Example Usage
```julia
𝕟 = lognormal([[100, 100, 1.1]]; d1 = 10.0, d2 = 1000.0, bins = 256)  
μ,σ = 80.0, 2000.0
T = (1.0 .+ erf.((𝕟.Dp .- μ)./(sqrt(2σ)))) 
𝕩 = T ⋅ 𝕟  
```
"""
function LinearAlgebra.:⋅(A::Vector{<:AbstractFloat}, 𝕟::SizeDistribution)
    if 𝕟.Dp[1] > 𝕟.Dp[2]
        nDp = reverse(A .* 𝕟.Dp)
        itpN = interpolate((nDp,), reverse(𝕟.N), Gridded(Linear()))
        extN = extrapolate(itpN, 0)

        itpS = interpolate((nDp,), reverse(𝕟.S), Gridded(Linear()))
        extS = extrapolate(itpS, 0)
        N = clean(extN(reverse(𝕟.Dp)))
        S = clean(extS(reverse(𝕟.Dp)))

        return SizeDistribution([[]], 𝕟.De, 𝕟.Dp, 𝕟.ΔlnD, reverse(S), reverse(N), :axdist)
    else
        nDp = A .* 𝕟.Dp
        itpN = interpolate((nDp,), reverse(𝕟.N), Gridded(Linear()))
        extN = extrapolate(itpN, 0)

        itpS = interpolate((nDp,), reverse(𝕟.S), Gridded(Linear()))
        extS = extrapolate(itpS, 0)
        N = clean(extN(𝕟.Dp))
        S = clean(extS(𝕟.Dp))

        return SizeDistribution([[]], 𝕟.De, 𝕟.Dp, 𝕟.ΔlnD, S, N, :axdist)
    end
end

⋅(𝕟::SizeDistribution, A::Vector{<:AbstractFloat}) =
    ⋅(A::Vector{<:AbstractFloat}, 𝕟::SizeDistribution)

"""
    +(𝕟₁::SizeDistribution, 𝕟₂::SizeDistribution)

Defines the sum of two size distributions. If diameter grids are not equal, then the
diameter grid of n2 is interpolated onto the n1 grid prior to addition.

```julia
𝕩 = 𝕟₁ + 𝕟₂ 
```
is defined such that 

```julia
𝕩.S = 𝕟₁.S + 𝕟₂.S 
𝕩.N = 𝕩.S .* 𝕟.ΔlnD 
```

Example Usage
```julia
𝕟₁ = lognormal([[120, 90, 1.20]]; d1 = 10.0, d2 = 1000.0, bins = 256)
𝕟₂ = lognormal([[90, 140, 1.15]]; d1 = 20.0, d2 = 800.0, bins = 64)
𝕩 = 𝕟₁ + 𝕟₂
```
"""
function +(𝕟₁::SizeDistribution, 𝕟₂::SizeDistribution)
    if 𝕟₁.Dp ≠ 𝕟₂.Dp
        itp = interpolate((𝕟₂.Dp,), 𝕟₂.N, Gridded(Linear()))
        ext = extrapolate(itp, 0)
        N = clean(ext(𝕟₁.Dp))

        itp = interpolate((𝕟₂.Dp,), 𝕟₂.S, Gridded(Linear()))
        ext = extrapolate(itp, 0)
        S = clean(ext(𝕟₁.Dp))
        S = 𝕟₁.S + S
    else
        S = 𝕟₁.S + 𝕟₂.S
    end
    N = S .* 𝕟₁.ΔlnD
    return SizeDistribution([[]], 𝕟₁.De, 𝕟₁.Dp, 𝕟₁.ΔlnD, S, N, :distsum)
end

"""
    -(𝕟₁::SizeDistribution, 𝕟₂::SizeDistribution)

Defines the sum of two size distributions. If diameter grids are not equal, then the
diameter grid of n2 is interpolated onto the n1 grid prior to addition.

```julia
𝕩 = 𝕟₁ - 𝕟₂ 
```
is defined such that 

```julia
𝕩.S = 𝕟₁.S - 𝕟₂.S 
𝕩.N = 𝕩.S .* 𝕟.ΔlnD 
```

Example Usage
```julia
𝕟₁ = lognormal([[120, 90, 1.20]]; d1 = 10.0, d2 = 1000.0, bins = 256)
𝕟₂ = lognormal([[90, 140, 1.15]]; d1 = 20.0, d2 = 800.0, bins = 64)
𝕩 = 𝕟₁ - 𝕟₂
```
"""
function -(𝕟₁::SizeDistribution, 𝕟₂::SizeDistribution)
    if 𝕟₁.Dp ≠ 𝕟₂.Dp
        itp = interpolate((𝕟₂.Dp,), 𝕟₂.N, Gridded(Linear()))
        ext = extrapolate(itp, 0)
        N = clean(ext(𝕟₁.Dp))

        itp = interpolate((𝕟₂.Dp,), 𝕟₂.S, Gridded(Linear()))
        ext = extrapolate(itp, 0)
        S = clean(ext(𝕟₁.Dp))
        S = 𝕟₁.S - S
    else
        S = 𝕟₁.S - 𝕟₂.S
    end
    N = S .* 𝕟₁.ΔlnD
    return SizeDistribution([[]], 𝕟₁.De, 𝕟₁.Dp, 𝕟₁.ΔlnD, S, N, :distsum)
end

"""
    interpolateDataFrameOntoδ(kw)

This function takes some measured size distribution in a DataFrame and interpolates 
it onto a DMA grid. kw is a tuple containing a DataFrame, symbols to columns to extract
which contain diameter and response function, and a DMA grid.

Example Usage
```julia
    𝕣 = (df, :Dp, :R, δ) |> interpolate_df_onto_thisδ
```

This extracts the columns Dp and R from df and interpolates it ont grid δ and
returns the results as a SizeDistribution. The df has to be sorted in ascending order.
"""
function interpolateDataFrameOntoδ(kw)
    df, δ  = kw[1], kw[end]
    Dp, Rcn = df[!,kw[2]], df[!,kw[3]]
    itp = interpolate((Dp,), Rcn, Gridded(Linear()))
    etp = extrapolate(itp, 0.0) 
    R = etp(δ.Dp)
    
    return SizeDistribution([],δ.De,δ.Dp,δ.ΔlnD,R./δ.ΔlnD,R,:interpolated)
end

"""
    interpolateSizeDistributionOntoδ(kw)

This function takes a size distribution and interpolates it onto a DMA grid. 
kw is a tuple containing a SizeDistribution and a DMA grid.

Example Usage
```julia
    𝕣 = (𝕟, δ) |> interpolateSizeDistributionOntoδ
```

This extracts the columns Dp and R from df and interpolates it ont grid δ and
returns the results as a SizeDistribution. The df has to be sorted in ascending order.
"""
function interpolateSizeDistributionOntoδ(kw)
    𝕟, δ = kw[1], kw[2]
    itp = interpolate((reverse(𝕟.Dp),), reverse(𝕟.S), Gridded(Linear()))
    etp = extrapolate(itp, 0.0) 
    S = etp(δ.Dp)
    
    return SizeDistribution([],δ.De,δ.Dp,δ.ΔlnD,S,S.*δ.ΔlnD,:interpolated)
end
