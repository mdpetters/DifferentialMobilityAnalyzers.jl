# Data Types

Three composite data types abstract the DMA setup. The type DMAconfig includes geometry, 
flow rate, and polarity.  The type DifferentialMobilityAnalyzer includes the DMA 
transmission functions, convolution matrices, and a native DMA mobility grid discretization. 
The type SizeDistribution includes a list of vectors to represent aerosol size distributions. 
Constructor functions are available to initialize these types.

## Index
```@index
Pages = ["types.md"]
```

## Types

### DMAconfig
```@docs
DMAconfig
```

### DifferentialMobilityAnalyzer
```@docs
DifferentialMobilityAnalyzer
```

### SizeDistribution
```@docs
SizeDistribution
```

### Regvars
```@docs
Regvars
```

## Constructor Functions

### setupDMA
```@docs
setupDMA
```

### setupSMPS
```@docs
setupSMPS
```

### setupSMPSdata
```@docs
setupSMPSdata
```

### lognormal
```@docs
lognormal
```

### DMALognormalDistribution
```@docs
DMALognormalDistribution    
```

### triangular
```@docs
triangular    
```

