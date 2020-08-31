# Physics

The following functions are used to define the physics of the instrument. 

## Index
```@index
Pages = ["physics.md"]
```

## Functions
### Fluid Viscosity
```@docs
η(Λ::DMAconfig)
```

### Cunningham Correction Factor
```@docs
cc(Λ::DMAconfig, d)
```

### Diffusion Coefficient
```@docs
dab(Λ::DMAconfig, d)
```

### dtoz
```@docs
dtoz(Λ::DMAconfig, d)
```

### ztod
```@docs
ztod
```

### vtoz
```@docs
vtoz(Λ::DMAconfig, v)
```

### ztov
```@docs
ztov(Λ::DMAconfig, z)
```

### Transfer Function
```@docs
Ω(Λ::DMAconfig, Z, zs)
```

```@docs
Ωav(Λ::DMAconfig, i::Int, k::Int; nint = 20)
```

### Charging Probability
```@docs
getTc(Λ::DMAconfig)
```

### Transmission Loss
```@docs
Tl(Λ::DMAconfig, Dp)
```