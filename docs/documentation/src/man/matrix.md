# Creating Convolution Matrices

By convention, bold upper case letters are used to denote all matrices. In code we use Unicode UTF-8 U+1D400 to U+1D419. It's easiest to copy paste the character from  a [UTF table](https://www.w3.org/TR/xml-entity-names/1D4.html). 

!!! note
    **UTF-8 Captial Bold Letters:** 𝐀, 𝐁, 𝐂, 𝐃, 𝐄, 𝐅, 𝐆, 𝐇, 𝐈, ...

## Fredholm Integral
The response in channel i (corresponding to a transmission through the DMA at a single voltage) is given by the Fredholm convolution integral 

 ``R_i = \int_{z_a}^{z_b} \sum_{k=1}^m \Omega(z,Z_{i,k}^s)T_c(D[z,1])T_l(D[z,1])\frac{dN}{d\ln D}\frac{d \ln D}{dz}dz \;\;\;\;\;\;\ i = 1,2...,n``
 
The integral is performed over the limits ``z_a`` and ``z_b``, which corresponds to the upper and lower mobility limit set by the voltage range used to operate the DMA. The function ``\frac{dN}{d\ln D}\frac{d\ln D}{dz}dz`` evaluates to the number concentration of particles that lie in the interval ``[z,z + dz]``. Note that ``D[z,1]`` is used in the charge filter and loss function since the integral is applied over the transform of the selected centroid mobility ``Z_{i,k}^s``. ``Z`` is the mobility vector of the grid and the subscript ``i`` denotes the response channel. The convolution integral can be discretized:  

 ``R_i = \sum_{j=1}^n \left [ \sum_{k=1}^m \Omega(Z_j,Z_{i,k}^s)T_c(D[Z_j,1])T_l(D[Z_j,1])N(Z_j) \right]``

``N(Z_j)`` is the the number concentration of particles in the ``j^{th}`` bin, ``i = 1...n``
are indices the observed instrument channel, ``j = 1...n`` are indices of the physical size bins, and ``k = 1...m`` are indices of charges carried by the particle. Here it is assumed that ``\Omega(Z_j,Z_{i,k}^s)`` can be approximated being constant over the bin ``[Ze_{j},Ze_{j+1}]``, which is only true if a sufficiently large number of size bins is used. 

The discretized version represents a set of ``n`` equations that can be written in matrix form 
    
``R = \mathbf{A}N``

where ``R`` is the response vector, ``\mathbf{A}`` is the convolution matrix, and ``N`` is the discretized true number concentration.

## Computing the Convolution Matrix
The terms comprising the convolution matrix 

``\Omega(Z_1,Z_{1,k}^s)T_c(D[Z_1,1])T_l(D[Z_1,1])``
    
can be written as 
    
```julia
T(zˢ,k,Λ,δ) = δ.Ω(Λ,δ.Z,zˢ/k).*δ.Tc(k,δ.Dp).*δ.Tl(Λ,δ.Dp)
```

The generic solution to constructing the convolution matrix is to define a forward transmission model and collect the terms. The resulting matrix is ``n \times n`` square matrix, where ``n`` equals to the number of bins.

```@example
using DifferentialMobilityAnalyzers              # hide
qsa,qsh = 1.66e-5, 8.33e-5                       # hide
t,p = 295.15, 1e5                                # hide
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369               # hide
leff = 13.0                                      # hide
m = 6                                            # hide
DMAtype = :cylindrical                           # hide
polarity = :-                                    # hide
Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,leff,polarity,m,DMAtype)  # hide
bins,z₁,z₂ = 60, vtoz(Λ,10000), vtoz(Λ,10)       
δ = setupDMA(Λ, z₁, z₂, bins)                   
T(zˢ,k,Λ,δ) = δ.Ω(Λ,δ.Z,zˢ/k).*δ.Tc(k,δ.Dp).*δ.Tl(Λ,δ.Dp)
𝐀 = (hcat(map(zˢ->Σ(k->T(zˢ,k,Λ,δ),Λ.m),δ.Z)...))'
```

See Notebook S2 in the [Notebooks](@ref) section for a step-by-step derivation. For a narrated description check out Session 2 of the [Tutorial](@ref). 

## Precomputed Matrices

The following matrices are precomputed for each DMA and stored in δ. The structure is mutable and the matrices can be altered if needed.

### Matrix 𝐀

The matrix 𝐀 describes transmission through the DMA with neutralizer and transmission loss function. The matrix 𝐀 is useful for modeling the measured response function of a size distribution passing through a stepping or scanning DMA.

```julia
δ = setupDMA(Λ, z₁, z₂, bins)                   
T(zˢ,k,Λ,δ) = δ.Ω(Λ,δ.Z,zˢ/k).*δ.Tc(k,δ.Dp).*δ.Tl(Λ,δ.Dp)
𝐀 = (hcat(map(zˢ->Σ(k->T(zˢ,k,Λ,δ),Λ.m),δ.Z)...))'
```

The type δ is mutable. You can therefore override the precomputed matrix by creating your own and manually adding it to the DMA grid.

```julia
T(zˢ,k,Λ,δ) = δ.Ω(Λ,δ.Z,zˢ/k).*δ.Tc(k,δ.Dp)
my𝐀 = (hcat(map(zˢ->Σ(k->T(zˢ,k,Λ,δ),Λ.m),δ.Z)...))'
δ.𝐀 .= my𝐀
```

### Matrix 𝐎

The matrix 𝐎 describes transmission through the DMA without neutralizer and with a transmission loss function. The matrix 𝐎 is useful for modeling the measured response function of a known mobility distribution passing through the second DMA in tandem DMA setups.

```julia
T(zˢ,k,Λ,δ) = δ.Ω(Λ,δ.Z,zˢ/k).*δ.Tl(Λ,δ.Dp)
𝐎 = (hcat(map(zˢ->Σ(k->T(zˢ,k,Λ,δ),Λ.m),δ.Z)...))'
```

### Matrix 𝐒

Talukdar and Swihart (2003) introduced the matrix 𝐒: sum the rows of 𝐀 and place the results on the diagonal of 𝐒. The matrix 𝐒 is used to compute an initial guess to constrain the Tikhonov inverse. 

```julia
𝐒 = zeros(bins, bins)
for i = 1:bins
	@inbounds 𝐒[i, i] = sum(𝐀[i, :])
end
```

### Matrix 𝐈

The identity matrix is used as weights matrix when computing the Tikhoniv inverse.

```julia
𝐈 = Matrix{Float64}(I, bins, bins)
```