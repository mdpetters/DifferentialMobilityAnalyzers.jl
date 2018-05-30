using Plots, DifferentialMobilityAnalyzers
plotlyjs()
T,p = 293.15, 1e5
qsa,β = 1.6666666e-5, 1/4
r1,r2,l = 9.37e-3,1.961e-2,0.44369
leff = 13.0
n = 6.0
Λ = DMAconfig(T,p,qsa,qsa/β,r1,r2,l,leff,:-,n)

v1,v2 = 10,10000
tscan, tc = 120,1
bins = tscan
z1,z2 = vtoz(Λ,v2), vtoz(Λ,v1)
δ = setupDMA(Λ, z1, z2, bins);
δ1 = setupSMPS(Λ, v1, v2, tscan, tc);
zˢ = dtoz(Λ, 200e-9)
i = 2
zˢ = δ.Z[i]
k = 2
plot(δ.Dp, δ.Ω(Λ,δ.Z,zˢ/k), xaxis = :log10)
plot!(δ1.Dp, δ1.Ω(Λ,i,k))
println(sum(δ.Ω(Λ,δ.Z,zˢ/k)))
println(sum(δ1.Ω(Λ,i,k)))
𝕟 = DMALognormalDistribution([[400, 30, 1.2],[500, 110, 1.7]], δ)
𝕟1 = DMALognormalDistribution([[400, 30, 1.2],[500, 110, 1.7]], δ1)
plot(δ1.Dp, δ1.𝐎*𝕟1.N, xscale = :log10)
plot!(δ.Dp, δ.𝐎*𝕟.N)
gui()
