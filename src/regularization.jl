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

# Function to setup the problem
function setupRegularization(𝐀, 𝐈, B, X₀)
    global Ψ = Regvars(𝐀, 𝐈, B, X₀)
end

# This function returns the inverted distribution as well as the
# L1 and L2 norms to construct the L-curve. The type of return
# value is optional, to facilitate definition of derivatives
function reginv(λs; r = :L1)
    Nλ = Array{Array{Float64}}(undef, 0)
    L1, L2 = Float64[], Float64[]
    for λ in λs
        Nx = inv(Ψ.𝐀' * Ψ.𝐀 + λ^2.0 * Ψ.𝐈) * (Ψ.𝐀' * Ψ.B + λ^2.0 * Ψ.X₀)
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

# Define the functions η, ρ and their derivatives. The functions
# are used to compute the curvature of the L-curve as defined in
# Eq.(14) of Hansen (2000)
η⁰(λ) = (log.(reginv(λ; r = :L2) .^ 2))[1]
ρ⁰(λ) = (log.(reginv(λ; r = :L1) .^ 2))[1]
ηᵖ(λ) = (derivative(η⁰, λ))[1]
ρᵖ(λ) = (derivative(ρ⁰, λ))[1]
ρ²ᵖ(λ) = (second_derivative(ρ⁰, λ))[1]
η²ᵖ(λ) = (second_derivative(η⁰, λ))[1]
κ(λ) = 2.0 * (ρᵖ(λ) * η²ᵖ(λ) - ηᵖ(λ) * ρ²ᵖ(λ)) / (ρᵖ(λ)^2.0 + ηᵖ(λ)^2.0)^1.5

# Compute the L-curve for n points between limits λ₁ and λ₂
function lcurve(λ₁::AbstractFloat, λ₂::AbstractFloat; n::Int = 10)
    λs = 10.0 .^ range(log10(λ₁), stop = log10(λ₂), length = n)
    L1, L2 = reginv(λs, r = :L1L2)
    κs = map(λ -> κ(λ), λs)
    ii = argmax(κs)
    if ii == length(κs)
        ii = ii - 1
    elseif ii == 1
        ii = ii + 1
    end
    return L1, L2, λs[ii-1:ii+1], ii
end

# Find the corner of the L-curve using iterative adjustment of the grid
function lcorner(λ₁::AbstractFloat, λ₂::AbstractFloat; n::Int = 10, r::Int = 3)
    L1, L2, λs, ii = lcurve(λ₁, λ₂; n = 10)
    for i = 1:r
        L1, L2, λs, ii = lcurve(λs[1], λs[3]; n = 10)
    end
    return λs[2]
end

# Warpper for the regularized inversion
function rinv(R, δ; λ₁ = 1e-2, λ₂ = 1e1)
    setupRegularization(δ.𝐀, δ.𝐈, R, inv(δ.𝐒) * R)  # setup the system
    λopt = lcorner(λ₁, λ₂; n = 10, r = 3)           # compute the optimal λ
    N = clean((reginv(λopt, r = :Nλ))[1])           # find the inverted size
    return SizeDistribution([], δ.De, δ.Dp, δ.ΔlnD, N ./ δ.ΔlnD, N, :regularized)
end
