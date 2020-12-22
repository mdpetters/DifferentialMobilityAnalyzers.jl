using RegularizationTools
using Underscores

t, p = 295.15, 1e5
qsa, qsh = 1.66e-5, 8.33e-5
r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369
Λ₁ = DMAconfig(t, p, qsa, qsh, r₁, r₂, l, 0.0, :-, 6, :cylindrical)
Λ₂ = DMAconfig(t, p, qsa, qsh, r₁, r₂, l, 0.0, :-, 6, :cylindrical)
𝕟 = lognormal([[9., 40., 1.5], [500., 180., 1.4]]; d1 = 800, d2 = 10.0, bins = 60)

Dd = 100e-9
dma2range = (Dd, 0.8, 3.0, 100)
mgf = 0.8:0.05:2.5 |> collect

f = TDMA1Ddomainfunction(𝕟, Λ₁, Λ₂, dma2range)
PDF(i) = @_ map(_ == i ? 1.0 : 0.0, 1:length(mgf))
A = designmatrix(mgf, f)
R1 = forwardmodel(mgf, PDF(10), mgf, f)
R2 = A*(PDF(10))
model = TDMA1Dpdf(𝕟, Λ₁, Λ₂, dma2range)
sd = model(𝕟, PDF(10), Dd, mgf)

@test R1 == R2
@test round(maximum(sd.N), digits = 0) == round(maximum(R2), digits = 0)

zˢ = dtoz(Λ₁, 100e-9) 
truegf = 1.6          
@test round(gfₖ(Λ₁, zˢ, truegf, 1), digits = 2) == 1.6
@test round(gfₖ(Λ₁, zˢ, truegf, 2), digits = 2) == 1.54
@test round(gfₖ(Λ₁, zˢ, truegf, 3), digits = 2) == 1.51
@test round(gfₖ(Λ₁, zˢ, truegf, 4), digits = 2) == 1.48  