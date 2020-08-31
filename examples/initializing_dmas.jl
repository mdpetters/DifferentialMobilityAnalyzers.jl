# Initializing DMAs


using DifferentialMobilityAnalyzers              

qsa,qsh = 1.66e-5, 8.33e-5                       # Qsample [m3 s-1], Qsheath [m3 s-1]
t,p = 295.15, 1e5                                # Temperature [K], Pressure [Pa]
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369               # DMA geometry [m]
leff = 13.0                                      # DMA effective diffusion length [m]
m = 6                                            # Upper number of charges to consider
DMAtype = :cylindrical                           # specify DMA type as cylindrical or radial
polarity = :-                                    # negative :- or positive :+ polartiy

Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,leff,polarity,m,DMAtype)  
bins,z₁,z₂ = 60, vtoz(Λ,10000), vtoz(Λ,10)       # bins, upper, lower mobility limit
δ = setupDMA(Λ, z₁, z₂, bins);                   # Setup DMA grid
