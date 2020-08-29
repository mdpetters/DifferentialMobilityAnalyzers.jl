# Operators

Operators are used to transform size distributions. The operators overload the Base or 
LinearAlgebra methods. Operators fall into two broad categories: operators changing 
number concentration and spectral density fields (𝕟.N and 𝕟.S) and operators that change 
the sizing vector (𝕟.Dp). The former include 𝕟₁ + 𝕟₂, 𝕟₁ - 𝕟₂, 𝕟₁ ∗ 𝕟₂, 𝕟₁ / 𝕟₂, a ∗ 𝕟, 
T .∗ 𝕟, and A ∗ 𝕟, while the latter include a · 𝕟 and T · 𝕟.

|  Operator |  Description |
|---|---|
| 𝕟₁ + 𝕟₂ | Superposition of the distributions 𝕟₁ and 𝕟₂ |
| 𝕟₁ - 𝕟₂ | Superposition of the distributions 𝕟₁ and 𝕟₂ |
| a ∗ 𝕟 | Uniform scaling of the concentration fields by factor a |
| 𝐀 ∗ 𝕟 | Matrix multiplication of 𝐀 and concentration |
| 𝕟₁ * 𝕟₂ | Scaled such that total number concentration equals to N1 ∗ N2 |
| 𝕟₁ / 𝕟₂ | Ratio of concentration fields of distributions 𝕟₁.N and 𝕟₂.N |
| a · 𝕟 | Uniform scaling of the diameter field of the size distribution by factor a |
| T · 𝕟 | Elementwise scaling the diameter field by factor T |
|||

## Index
```@index
Pages = ["operators.md"]
```

## Number Operators

```@docs
+(𝕟₁::SizeDistribution, 𝕟₂::SizeDistribution)
```

```@docs
-(𝕟₁::SizeDistribution, 𝕟₂::SizeDistribution)
```

```@docs
*(a::Number, 𝕟::SizeDistribution)
```

```@docs
*(𝐀::AbstractMatrix, 𝕟::SizeDistribution)
```

```@docs
*(𝕟₁::SizeDistribution, 𝕟₂::SizeDistribution)
```

```@docs
/(𝕟₁::SizeDistribution, 𝕟₂::SizeDistribution)
```



## Size Operators

```@docs
LinearAlgebra.:⋅(a::Number, 𝕟::SizeDistribution)
```

```@docs
LinearAlgebra.:⋅(A::Vector{<:AbstractFloat}, 𝕟::SizeDistribution)
```