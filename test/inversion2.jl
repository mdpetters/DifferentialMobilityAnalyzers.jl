using Random

# Test inversion from Notebook 4
t, p = 295.15, 1e5
qsa, β = 1.66e-5, 1 / 5
r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369
leff = 13.0
m = 3
Λ = DMAconfig(t, p, qsa, qsa / β, r₁, r₂, l, leff, :-, m, :cylindrical)
bins, z₁, z₂ = 128, dtoz(Λ, 1000e-9), dtoz(Λ, 10e-9)
δ = setupDMA(Λ, z₁, z₂, bins);

𝐀 =
    (hcat(map(
        i -> Σ(k -> δ.Ω(Λ, δ.Z, i / k) .* δ.Tc(k, δ.Dp) .* δ.Tl(Λ, δ.Dp), Λ.m),
        δ.Z,
    )...))'

𝕟 = DMALognormalDistribution([[400, 30, 1.2], [500, 110, 1.7]], δ)
Random.seed!(703)
tscan = 120
Qcpc = 16.66
t = tscan ./ bins

𝕣 = δ.𝐀 * 𝕟;

𝕣1 = 𝐀 * 𝕟;
@test 𝕣.N == 𝕣1.N

c = 𝕣.N * Qcpc * t;
R = Float64[]
for i in c
    f = rand(Poisson(i), 1)
    push!(R, f[1] / (Qcpc * t))
end

λ₁, λ₂ = 1e-3, 1e1
eyeM = Matrix{Float64}(I, bins, bins)
setupRegularization(δ.𝐀, eyeM, R, inv(δ.𝐒) * R, 1)
λopt = lcorner(λ₁, λ₂; n = 10, r = 3)
N = clean((reginv(λopt, r = :Nλ))[1])
𝕟ᵢₙᵥ = SizeDistribution([], 𝕟.De, 𝕟.Dp, 𝕟.ΔlnD, N ./ 𝕟.ΔlnD, N, :regularized)
@test round(Int, sum(𝕟ᵢₙᵥ.N)) == 890
