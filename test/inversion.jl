# Test inversion from Notebook 4
t,p = 295.15, 1e5
qsa,β = 1.66e-5, 1/5
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369
leff = 13.0
m = 3
Λ = DMAconfig(t,p,qsa,qsa/β,r₁,r₂,l,leff,:-,m)
bins,z₁,z₂ = 128, dtoz(Λ,1000e-9), dtoz(Λ,10e-9)
δ = setupDMA(Λ, z₁, z₂, bins);
𝕟 = DMALognormalDistribution([[400, 30, 1.2],[500, 110, 1.7]], δ)
srand(703)
tscan = 120
Qcpc = 16.66
t = tscan./bins

𝐀=(hcat(map(i->Σ(k->δ.Ω(Λ,δ.Z,i/k).*δ.Tc(k,δ.Dp).*δ.Tl(Λ,δ.Dp),Λ.m),δ.Z)...))'
@test round.(sum(𝐀),2) == 77.38
@test 𝐀 == δ.𝐀

T = (zˢ,k,Λ,δ) -> δ.Ω(Λ,δ.Z,zˢ/k).*δ.Tc(k,δ.Dp).*δ.Tl(Λ,δ.Dp)
𝐀=(hcat(map(zˢ->Σ(k->T(zˢ,k,Λ,δ),Λ.m),δ.Z)...))'
@test round.(sum(𝐀),2) == 77.38

𝕣 = δ.𝐀*𝕟;
c = 𝕣.N*Qcpc*t;
R = Float64[]
for i = c
    f = rand(Poisson(i),1)
    push!(R,f[1]/(Qcpc*t))
end
λ₁,λ₂ = 1e-3, 1e1
setupRegularization(δ.𝐀,eye(bins),R,inv(δ.𝐒)*R)
λopt = lcorner(λ₁,λ₂;n=10,r=3)
N =  clean((reginv(λopt, r = :Nλ))[1])
𝕟ᵢₙᵥ= SizeDistribution([],𝕟.De,𝕟.Dp,𝕟.ΔlnD,N./𝕟.ΔlnD,N,:regularized)
@test round(Int, sum(𝕟ᵢₙᵥ.N)) == 894
