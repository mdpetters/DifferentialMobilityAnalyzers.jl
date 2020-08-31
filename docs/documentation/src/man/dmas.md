# Initializing DMAs

## Notation
This package is designed to work with chained DMA systems. To isolate individual DMA systems, the instruments are abstracted into the data types [DMAconfig](@ref) and [DifferentialMobilityAnalyzer](@ref). By convention variables of type DMAconfig are assigned the letter Λ and variables of type DifferentialMobilityAnalyzer the letter δ. Subscripts or superscripts can be used to distinguish DMA1, DMA2, ... That is Λ₁, Λ₂, ... and δ₁, δ₂, ...

## DMA Configuration
The DMA is an annular capacitor. The column's properties are defined by the radii $r_1$, $r_2$, the length of the aerosol path, $l$. Operation conditions are defined by the the electric potential $v$ applied across the annulus and the four flow rates: sheath flow, $q_{sh}$, polydisperse aerosol flow $q_a$, excess flow, $q_{ex}$, and monodisperse sample flow, $q_{sa}$. Throughout this work it is assumed that the flows are balanced, i.e. $q_{sh} = q_{ex}$ and $q_{sa} = q_a$. The two flows tracked are $q_{sh}$ and $q_{sa}$.

Here is an example initialization of the DMAconfig

```@example
using DifferentialMobilityAnalyzers              

qsa,qsh = 1.66e-5, 8.33e-5                       # Qsample [m3 s-1], Qsheath [m3 s-1]
t,p = 295.15, 1e5                                # Temperature [K], Pressure [Pa]
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369               # DMA geometry [m]
leff = 13.0                                      # DMA effective diffusion length [m]
m = 6                                            # Upper number of charges to consider
DMAtype = :cylindrical                           # specify DMA type as cylindrical or radial
polarity = :-                                    # negative :- or positive :+ polartiy

Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,leff,polarity,m,DMAtype)  
```

The data type Λ defines the DMA in terms of flow rates, geometry and power supply polarity. The geometry parameters used in this example correspond to the TSI 3080 long column. The effective diffusion length describes [Transmission Loss](@ref) in the DMA column. Diffusional loss is ignored if leff = 0. The DMA configuration is used to cross reference to the appropriate [Fluid Viscosity](@ref), [Cunningham Correction Factor](@ref), [Diffusion Coefficient](@ref), particle [Charging Probability](@ref), and DMA [Transfer Function](@ref) for the given thermodynamic state and fluid velocity. Check out the [Tutorial](@ref) Session 1 and/or Notebook S1 in the [Notebooks](@ref) section for visualizations.

## DMA Grid
The DMA grid encodes the actual operation of the instrument. The DMA is operates between an upper and lower voltage limit. The full range is usually 10V to 10kV. At higher voltages the electric field breaks down. A convenient way to bin the DMA is to work with a log spaced mobility grid, which in essence is the size distribution histogram where diameters are computed using [ztod](@ref). 

The DMA grid is instantiated using one of the contructor functions: [setupDMA](@ref), [setupSMPS](@ref), [setupSMPSdata](@ref). Each of these defines the lower and upper size limit and the number of bins of the grid. 

```@example
using DifferentialMobilityAnalyzers              # hide
qsa,qsh = 1.66e-5, 8.33e-5                       # hide
t,p = 295.15, 1e5                                # hide
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369               # hide
leff = 13.0                                      # hide
m = 6                                            # hide
DMAtype = :cylindrical                           # hide
polarity = :-                                    # hide
Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,leff,polarity,m,DMAtype) # hide  
bins,z₁,z₂ = 60, vtoz(Λ,10000), vtoz(Λ,10)    # bins, upper, lower mobility limit
δ = setupDMA(Λ, z₁, z₂, bins);                # Setup DMA grid
```

The resulting type δ contains 
- ```Z``` are the mobility bin midpoints, 
- ```Dp``` are the diameter bin midpoints (internally in units of nm) 
- ```Ze``` are the mobility bin edges, 
- ```De``` are the diameter bin edges, 
- ```ΔlnD``` is the natural log ratio of upper and lower diameter bin edge
- the matrices 𝐀, 𝐒, 𝐎, 𝐈 used in the size distribution [Inversion Routines](@ref)
- the function Ω ([Transfer Function](@ref)), Tc ([Charging Probability](@ref)), Tl ([Transmission Loss](@ref))

Check out the [Tutorial](@ref) Session 1 and/or Notebook S1 in the [Notebooks](@ref) section for visualizations.