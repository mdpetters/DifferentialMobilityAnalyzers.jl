# +
# This file handles the regularization for size distribution inversion
#
# Author: Markus Petters (mdpetter@ncsu.edu)
# 	      Department of Marine Earth and Atmospheric Sciences
#         NC State University
#         Raleigh, NC 27605
#
#         April, 2018
#-

@doc raw"""
    setupRegularization(𝐀, 𝐈, B, X₀, n)

Initialize the [Regvars](@ref) used to compute the Tikhonov regularization. Regvars 
is a data type that stores the inversion problem setup. It also stores the precomputed 
A'A matrix for performance optimization
- A is the convolution matrix
- I is the identity matrix
- B is the response vector
- X0 is the initial guess
- n is the number of BLAS threads 
"""
function setupRegularization(𝐀, 𝐈, B, X₀, n)
    global Ψ = Regvars(𝐀[:,:], 𝐈, B, X₀,(𝐀'𝐀)[:,:],n)
end

# This function returns the inverted distribution as well as the
# L1 and L2 norms to construct the L-curve. The type of return
# value is optional, to facilitate definition of derivatives

# Keep for backward compatibility
function reginv(λs; r = :L1)
    Nλ = Array{Array{Float64}}(undef, 0)
    L1, L2 = Float64[], Float64[]
    for λ in λs
        Nx = Ninv(λ)
        push!(L1, norm(Ψ.𝐀 * Nx - Ψ.B))
        push!(L2, norm(Ψ.𝐈 * (Nx - Ψ.X₀)))
        push!(Nλ, Nx)
    end
    if r == :L1
        return L1
    elseif r == :L2
        return L2
    elseif r == :L1L2
        return L1, L2
    elseif r == :Nλ
        return Nλ
    end
end

function zot(A::AbstractMatrix, λ::AbstractFloat) 
    a = deepcopy(A)
    n = size(a, 1)
    for i = 1:n
        @inbounds a[i,i] += λ
    end
    return a
end


@doc raw"""
    Ninv(λ)

This function computes the regularized inverse for a specified regularization parameter λ.
This requires that the problem is iniatialized throug [setupRegularization](@ref). 
Use the [clean](@ref) function to truncate negative values.

Example Usage
```julia
# R is the response vector
setupRegularization(δ.𝐀, δ.𝐈, R, inv(δ.𝐒) * R, n)  
N = clean(Ninv(0.5))                           
```
"""
Ninv(λ::AbstractFloat) = 
    cholesky!(Hermitian(zot(Ψ.AA, λ^2.0))) \ (Ψ.𝐀' * Ψ.B + λ^2.0 * Ψ.X₀)

@doc raw"""
    L1L2(λ::AbstractFloat)

Returns the L1 and L2 norm for regularization parameter λ

Example Usage
```julia
# R is the response vector
setupRegularization(δ.𝐀, δ.𝐈, R, inv(δ.𝐒) * R, n)  
L1, L2 = L1L2(0.5)
```
"""
function L1L2(λ::AbstractFloat)
    Nx = Ninv(λ) 
    return norm(Ψ.𝐀 * Nx - Ψ.B), norm(Ψ.𝐈 * (Nx - Ψ.X₀))
end

@doc raw"""
    L1(λ::AbstractFloat)

Returns the L1 norm for regularization parameter λ

Example Usage
```julia
# R is the response vector
setupRegularization(δ.𝐀, δ.𝐈, R, inv(δ.𝐒) * R, n)  
L1 = L1(0.5)
```
"""
function L1(λ::AbstractFloat)
    Nx = Ninv(λ)
    return norm(Ψ.𝐀 * Nx - Ψ.B)
end

@doc raw"""
    L2(λ::Float64)

Returns the L2 norm for regularization parameter λ

Example Usage
```julia
# R is the response vector
setupRegularization(δ.𝐀, δ.𝐈, R, inv(δ.𝐒) * R, n)  
L2 = L2(0.5)
```
"""
function L2(λ::AbstractFloat)
    Nx = Ninv(λ)
    return norm(Ψ.𝐈 * (Nx - Ψ.X₀))
end


# Define the functions η, ρ and their derivatives. The functions
# are used to compute the curvature of the L-curve as defined in
# Eq.(14) of Hansen (2000)
η⁰(λ) = (log.(L2.(λ) .^ 2))[1]
ρ⁰(λ) = (log.(L1.(λ) .^ 2))[1]
ηᵖ(λ) = (derivative(η⁰, λ))[1]
ρᵖ(λ) = (derivative(ρ⁰, λ))[1]
ρ²ᵖ(λ) = (second_derivative(ρ⁰, λ))[1]
η²ᵖ(λ) = (second_derivative(η⁰, λ))[1]
function κ(λ::AbstractFloat)  
    rᵖ = ρᵖ(λ)
    nᵖ = ηᵖ(λ)
    2.0 * (rᵖ * η²ᵖ(λ) - nᵖ * ρ²ᵖ(λ)) / (rᵖ^2.0 + nᵖ^2.0)^1.5
end

# Compute the L-curve for n points between limits λ₁ and λ₂
function lcurve(λ₁::AbstractFloat, λ₂::AbstractFloat; n::Int = 10)
	BLAS.set_num_threads(Ψ.n)
    λs = 10.0 .^ range(log10(λ₁), stop = log10(λ₂), length = n)
    L1, L2 = L1L2.(λs)
    κs = map(λ -> κ(λ), λs)
    ii = argmax(κs)
    if ii == length(κs)
        ii = ii - 1
    elseif ii == 1
        ii = ii + 1
    end
    return L1, L2, λs[ii-1:ii+1], ii
end

@doc raw"""
    lcorner(λ₁::AbstractFloat, λ₂::AbstractFloat; n::Int = 10, r::Int = 3)

Searches for the corner of the L-curve between λ₁ and λ₂. 
- n is the number of λs to divide between between λ₁ and λ₂. The alogorithm picks the best 
λ from these ten. Then it repeats for a narrower range of λ₁ and λ₂ around this best value
- r is the number of times to repeat the search

Example Usage
```julia
# R is the response vector
setupRegularization(δ.𝐀, δ.𝐈, R, inv(δ.𝐒) * R, n)  # setup the system
λopt = lcorner(λ₁, λ₂; n = 10, r = 3)              # compute the optimal λ
```
"""
function lcorner(λ₁::AbstractFloat, λ₂::AbstractFloat; n::Int = 10, r::Int = 3)
    L1, L2, λs, ii = lcurve(λ₁, λ₂; n = n)
    for i = 1:r
        L1, L2, λs, ii = lcurve(λs[1], λs[3]; n = n)
    end
    return λs[2]
end

@doc raw"""
    rinv(R::AbstractVector, δ::DifferentialMobilityAnalyzer; λ₁ = 1e-2, λ₂ = 1e1, n = 1)

The function rinv is a wrapper to perform the Tikhonov inversion.
- R the response vector to be inverted
- δ the DifferentialMobilityAnalyzer from which R was produced
- λ₁ and λ₂ are the bounds of the search for the optical regularization parameter
- n is the number of BLAS threads to use (currently 1 is fastest)

The function returns an inverted size distribution of type [SizeDistribution](@ref)

Example Usage 

```julia
# Load Data
df = CSV.read("example_data.csv", DataFrame)

# Setup the DMA
t, p, lpm = 293.15, 940e2, 1.666e-5      
r₁, r₂, l = 9.37e-3,1.961e-2,0.44369     
Λ = DMAconfig(t,p,1lpm,4lpm,r₁,r₂,l,0.0,:+,6,:cylindrical)  
δ = setupDMA(Λ, vtoz(Λ,10000), vtoz(Λ,10), 120)

# Interpolate the data onto the DMA grid
𝕣 = (df, :Dp, :Rcn, δ) |> interpolateDataFrameOntoδ

# Compute the inverse. 
𝕟ⁱⁿᵛ = rinv(𝕣.N, δ, λ₁ = 0.1, λ₂ = 1.0)
```
"""
function rinv(R::AbstractVector, δ::DifferentialMobilityAnalyzer; λ₁ = 1e-2, λ₂ = 1e1, n = 1)
    setupRegularization(δ.𝐀, δ.𝐈, R, inv(δ.𝐒) * R, n)  # setup the system
    λopt = lcorner(λ₁, λ₂; n = 10, r = 3)           # compute the optimal λ
    N = clean(Ninv(λopt))                           # find the inverted size
    return SizeDistribution([], δ.De, δ.Dp, δ.ΔlnD, N ./ δ.ΔlnD, N, :regularized)
end
