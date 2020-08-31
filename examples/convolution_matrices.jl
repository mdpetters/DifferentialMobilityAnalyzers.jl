
using DifferentialMobilityAnalyzers

qsa, qsh = 1.66e-5, 8.33e-5                       # Qsample [m3 s-1], Qsheath [m3 s-1]
t, p = 295.15, 1e5                                # Temperature [K], Pressure [Pa]
r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369               # DMA geometry [m]
leff = 13.0                                      # DMA effective diffusion length [m]
m = 6                                            # Upper number of charges to consider
DMAtype = :cylindrical                           # specify DMA type as cylindrical or radial
polarity = :-                                    # negative :- or positive :+ polartiy

Λ = DMAconfig(t, p, qsa, qsh, r₁, r₂, l, leff, polarity, m, DMAtype)
bins, z₁, z₂ = 60, vtoz(Λ, 10000), vtoz(Λ, 10)       # bins, upper, lower mobility limit
δ = setupDMA(Λ, z₁, z₂, bins);                   # Setup DMA grid

# Convolution matrix A
T(zˢ, k, Λ, δ) = δ.Ω(Λ, δ.Z, zˢ / k) .* δ.Tc(k, δ.Dp) .* δ.Tl(Λ, δ.Dp)
𝐀 = (hcat(map(zˢ -> Σ(k -> T(zˢ, k, Λ, δ), Λ.m), δ.Z)...))'

# Convolution matrix O (same as A but without charge filter)
T(zˢ, k, Λ, δ) = δ.Ω(Λ, δ.Z, zˢ / k) .* δ.Tl(Λ, δ.Dp)
𝐎 = (hcat(map(zˢ -> Σ(k -> T(zˢ, k, Λ, δ), Λ.m), δ.Z)...))'

# Identity matrix 
𝐈 = Matrix{Float64}(I, bins, bins)

# Initial guess matrix
𝐒 = zeros(bins, bins)
for i = 1:bins
    @inbounds 𝐒[i, i] = sum(𝐀[i, :])
end

# Update precomputed matrix with custom matrix
T(zˢ, k, Λ, δ) = δ.Ω(Λ, δ.Z, zˢ / k) .* δ.Tc(k, δ.Dp)
my𝐀 = (hcat(map(zˢ -> Σ(k -> T(zˢ, k, Λ, δ), Λ.m), δ.Z)...))'
δ.𝐀 .= my𝐀
