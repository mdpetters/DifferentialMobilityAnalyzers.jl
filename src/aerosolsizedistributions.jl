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

# λ functions
md = (A, x) -> @. A[1] / (√(2π) * log(A[3])) * exp(-(log(x / A[2]))^2 / (2log(A[3])^2))
logn = (A, x) -> mapreduce((A) -> md(A, x), +, A)
clean = x -> map(x -> x < 0.0 ? 0.0 : x, x)

# Multimodal lognormal size distribution on a generic logspace grid
function lognormal(A; d1 = 8.0, d2 = 2000.0, bins = 256)
    #De = logspace(log10(d1), log10(d2), bins+1)
    De = 10.0 .^ range(log10(d1), stop = log10(d2), length = bins + 1)
    Dp = sqrt.(De[2:end] .* De[1:end-1])
    ΔlnD = log.(De[2:end] ./ De[1:end-1])
    S = logn(A, Dp)
    N = S .* ΔlnD
    return SizeDistribution(A, De, Dp, ΔlnD, S, N, :lognormal)
end

# Triangular size distribution for a given diameter and number concentration
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

# Lognormal distrubution on a DMA size grid
function DMALognormalDistribution(A, δ::DifferentialMobilityAnalyzer)
    S = logn(A, δ.Dp)

    return SizeDistribution(A, δ.De, δ.Dp, δ.ΔlnD, S, S .* δ.ΔlnD, :DMA)
end

# -------------------------- Size Distribution Arithmetic ---------------------
# --------------------------- Block 1: * .* and ./  ---------------------------
function *(a::AbstractFloat, 𝕟::SizeDistribution)
    # This function defines the product of a scalar and a size distribution
    N = a * 𝕟.N
    S = a * 𝕟.S

    return SizeDistribution([[]], 𝕟.De, 𝕟.Dp, 𝕟.ΔlnD, S, N, :axdist)
end

*(𝕟::SizeDistribution, a::AbstractFloat) = *(a::AbstractFloat, 𝕟::SizeDistribution)

function *(a::Vector{<:AbstractFloat}, 𝕟::SizeDistribution)
    # This function defines the product of a vector and a size distribution
    N = a .* 𝕟.N
    S = a .* 𝕟.S
    return SizeDistribution([[]], 𝕟.De, 𝕟.Dp, 𝕟.ΔlnD, S, N, :axdist)
end

*(𝕟::SizeDistribution, a::Vector{<:AbstractFloat}) =
    *(a::Vector{<:AbstractFloat}, 𝕟::SizeDistribution)


function *(𝐀::AbstractMatrix, 𝕟::SizeDistribution)
    # This function defines the product of a matrix and a size distribution
    N = 𝐀 * 𝕟.N
    S = 𝐀 * 𝕟.S
    return SizeDistribution([[]], 𝕟.De, 𝕟.Dp, 𝕟.ΔlnD, S, N, :axdist)
end

*(𝕟::SizeDistribution, 𝐀::AbstractMatrix) = *(𝐀::AbstractMatrix, 𝕟::SizeDistribution)

function *(𝕟₁::SizeDistribution, 𝕟₂::SizeDistribution)
    # This function defines the product of two size distributions
    Nsq = 𝕟₁.N .* 𝕟₂.N
    N = sum(𝕟₁.N) * sum(𝕟₂.N) * Nsq ./ sum(Nsq)
    S = N ./ 𝕟₁.ΔlnD
    return SizeDistribution([[]], 𝕟₁.De, 𝕟₁.Dp, 𝕟₁.ΔlnD, S, N, :dist_sq)
end


function /(𝕟₁::SizeDistribution, 𝕟₂::SizeDistribution)
    # This function defines the product of two size distributions
    N = 𝕟₁.N ./ 𝕟₂.N
    S = 𝕟₁.S ./ 𝕟₂.S
    return SizeDistribution([[]], 𝕟₁.De, 𝕟₁.Dp, 𝕟₁.ΔlnD, S, N, :dist_sq)
end


"""
    ⋅(a::AbstractFloat, 𝕟::SizeDistribution)

Multiplication of a scalare and a size distribution. The net result is a uniform diameter 
shift of the size distribution. The function is symmetric such that a⋅𝕟 == 𝕟⋅a

Let a denote a floating point scalar and 𝕟 denote a size distribution. Then
```julia
𝕩 = a⋅𝕟
```
is defined such that 
```julia
𝕩.Dp = a*𝕟.Dp 
```

Example Usage
```julia
a = 2.0 # Note that a must be a floating point number
𝕟 = lognormal([[300, 100, 1.3]]; d1 = 10.0, d2 = 1000.0, bins = 256)
𝕩 = a⋅𝕟 
```
"""
function LinearAlgebra.:⋅(a::AbstractFloat, 𝕟::SizeDistribution)
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

#function LinearAlgebra.:.⋅(A::Array{Float64,1}, 𝕟::SizeDistribution)
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
𝕟₂ = lognormal([[90, 140, 1.15]]; d1 = 20.0, d2 = 800.0, bins = 64
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

function -(𝕟₁::SizeDistribution, 𝕟₂::SizeDistribution)
    # This function defines the subtraction of two size distributions

    # If grids are not equal, then interpolate n2 onto n1 grid
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
    interpolate_df_onto_thisδ(kw)

This function takes some measured size distribution in a DataFrame and and interpolates 
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
    
    return SizeDistribution([],δ.De,δ.Dp,δ.ΔlnD,R./δ.ΔlnD,R,:response)
end