module DifferentialMobilityAnalyzers
# +
# This file defines the module DifferentialMobilityAnalyzers.jl
#
# Author: Markus Petters (mdpetter@ncsu.edu)
# 	      Department of Marine Earth and Atmospheric Sciences
#         NC State University
#         Raleigh, NC 27605
#
#         May, 2018
#-

using Interpolations
using SpecialFunctions
using DataFrames
using Calculus
using LinearAlgebra
using Random
using CSV
using Distributions
using RegularizationTools
using Lazy
using MLStyle
using Memoize
using Underscores
import Lazy.@>, Lazy.@>>, Lazy.@as

import Base.*,                # Import to extend operators
    Base./,                   # Import to extend operators
    Base.+,                   # Import to extend operators
    Base.-,                   # Import to extend operators
    LinearAlgebra.⋅

export DMAconfig,                    # Data type to hold DMA config
    DifferentialMobilityAnalyzer,    # Data type to hold grid and matrices
    SizeDistribution,         # Data type to hold size distribution
    Regvars,                  # Data type to hold regularization problem
    Ψ,                        # Data structure that contains regvars
    setupDMA,                 # Function to initialize DMA
    setupSMPS,                # Function to initialize SMPS
    setupSMPSdata,            # Function to initialize SMPS with data
    setupDMAgridded,          # Function to initialize SMPS with De
    setupRegularization,      # Function to initialize regularization
    clean,                    # Function to remove negative numbers
    Σ,                        # Function to sum functions
    λ,                        # Function to compute mean free path
    η,                        # Function to compute air viscosity
    Ωav,                      # Function to compute transfer through the DMA
    Ω,                        # Function to compute transfer through the DMA
    cc,                       # Function to compute slip-flow correction
    dab,                      # Function to compute diffusion coefficient
    vtoz,                     # Function to convert Voltage to Mobility
    ztov,                     # Function to convert Mobility to Voltage
    dtoz,                     # Function to convert Diameter to Mobility
    ztod,                     # Function to convert Mobility to Diameter
    lognormal,                # Function to compute lognormal size dist
    triangular,               # Function to compute DMA triangular size dist
    getTc,                    # Function to compute charge distribution function
    Tl,                       # Function to compute transmission efficiency
    DMALognormalDistribution, # Function to compute lognormal size dist
    reginv,                   # Function to compute regularized inverse
    Ninv,                     # Function to compute regularized inverse
    lcurve,                   # Function to compute the L-curve
    lcorner,                  # Function to find corner of L-curve
    rinv,                     # Function to compute the inverse
    rinv2,                    # Function to compute the inverse with RegularizationTools
    L1,                       # Function to compute L1 norm
    L2,                       # Function to compute L2 norm
    L1L2,                     # Function to compute L1 and L2 norms
    loadtsidata,              # Function to loads TSI data file; loadtsidata.
    β12brown,                 # Brownian coagualation kernel
    β12zebel,                 # Zebel charge enhancement
    β12,                      # Brownian+Zebel coagulation kernel
    figure,                   # Function to setup figures
    interpolateDataFrameOntoδ,# Interpolate data onto DMA grid
    interpolateSizeDistributionOntoδ,# Interpolate data onto DMA grid
    benchmark,                # Single benchmark ren
    runbenchmarks,            # Create summary of benchmarks
    gfₖ,                      # Effective growth factor multi-charge particles
    TDMA1Dpdf,              # 1D PDF model of TDMA transfer
    TDMA1Ddomainfunction,     # 1D domain function to compute design matrix    
    initializeDefaultMatrices, # Precompute inversion matrices
    Air,
    N2

abstract type CarrierGas end
struct Air <: CarrierGas end
struct N2 <: CarrierGas end

@doc raw"""
    DMAconfig

Data type to abstract the DMA geometry and state of the fluid. 

    t::AbstractFloat          # Temperature [K]
    p::AbstractFloat          # Pressure [Pa]
    qsa::AbstractFloat        # Sample flow [m3 s-1]
    qsh::AbstractFloat        # Sheath flow [m3 s-1]
    r1::AbstractFloat         # Inner column radius [m]
    r2::AbstractFloat         # Outer column radius [m]
    l::AbstractFloat          # Column length [m]
    leff::AbstractFloat       # Effective length [m]
    polarity::Symbol          # Power supply polarity [:+] or [:-]
    m::Int8                   # Number of charges in charge correction [-]
    DMAtype::Symbol           # Designate :radial, :cylindrical

Example Usage
```julia
t,p = 295.15, 1e5                        
qsa,qsh = 1.66e-5, 8.3e-5                     
r₁,r₂,l = 9.37e-3,1.961e-2,0.44369               
Λ = DMAconfig(t,p,qsa,qsh,r₁,r₂,l,13.0,:-,6,:cylindrical) 
```julia

!!! note

    When defining a radial DMA, r₁,r₂,l map to  r₁,r₂,b as defined in Zhang Shou-Hua Zhang, 
    Yoshiaki Akutsu, Lynn M. Russell, Richard C. Flagan & John H. Seinfeld (1995) 
    Radial Differential Mobility Analyzer, Aerosol Science and Technology, 
    23:3, 357-372, DOI: 10.1080/02786829508965320.

"""
struct DMAconfig
    t::AbstractFloat          # Temperature [K]
    p::AbstractFloat          # Pressure [Pa]
    qsa::AbstractFloat        # Sample flow [m3 s-1]
    qsh::AbstractFloat        # Sheath flow [m3 s-1]
    r1::AbstractFloat         # Inner column radius [m]
    r2::AbstractFloat         # Outer column radius [m]
    l::AbstractFloat          # Column length [m]
    leff::AbstractFloat       # Effective length [m]
    polarity::Symbol          # Power supply polarity [:+] or [:-]
    m::Int8                   # Number of charges in charge correction [-]
    DMAtype::Symbol           # Designate :radial, :cylindrical
    gas::CarrierGas           # CarrierGas, either Air() or N2()
end

DMAconfig(t, p, qsa, qsh, r1, r2, l, leff, polarity, m, DMAtype) = 
    DMAconfig(t, p, qsa, qsh, r1, r2, l, leff, polarity, m, DMAtype, Air())

@doc raw"""
    DifferentialMobilityAnalyzer

The type DifferentialMobilityAnalyzer contains the DMA transmission functions, 
a discretized mobility grid to represent the mobility distribution and precomputed
convolution matrices.

    Ω::Function                    # DMA transfer function
    Tc::Function                   # Charge filter Function
    Tl::Function                   # DMA Penetration efficiency function
    Z::Vector{<:AbstractFloat}     # Mobility array midpoints
    Ze::Vector{<:AbstractFloat}    # Mobility array bin edges
    Dp::Vector{<:AbstractFloat}    # Mobility diameter midpoints
    De::Vector{<:AbstractFloat}    # Mobility diameter bin edges
    ΔlnD::Vector{<:AbstractFloat}  # ln(de[i+1])-ln(de[i])
    𝐀::AbstractMatrix              # Convolution matrix
    𝐒::AbstractMatrix              # Convolution matrix for initial guess
    𝐎::AbstractMatrix              # Convolution matrix for no charge filter
    𝐈::AbstractMatrix               # IdentiyMatrix

The field is initialized using one of the the constructor functions:
- [setupDMA](@ref)
- [setupSMPS](@ref)
- [setupSMPSdata](@ref)

!!! note

    Diameters stored in δ are in units of nm. Mobility in m2 s-1 V-1. The functions 
    [Transfer Function](@ref) Ω, [Charging Probability](@ref) Tc and 
    [Transmission Loss](@ref) Tl take diameter in units of nm

"""
struct DifferentialMobilityAnalyzer
    Ω::Function                    # DMA transfer function
    Tc::Function                   # Charge filter Function
    Tl::Function                   # DMA Penetration efficiency function
    Z::Vector{<:AbstractFloat}     # Mobility array midpoints
    Ze::Vector{<:AbstractFloat}           # Mobility array bin edges
    Dp::Vector{<:AbstractFloat}           # Mobility diameter midpoints
    De::Vector{<:AbstractFloat}           # Mobility diameter bin edges
    ΔlnD::Vector{<:AbstractFloat}         # ln(de[i+1])-ln(de[i])
    𝐀::AbstractMatrix            # Convolution matrix
    𝐒::AbstractMatrix            # Convolution matrix for initial guess
    𝐎::AbstractMatrix            # Convolution matrix for no charge filter
    𝐈::AbstractMatrix            # IdentiyMatrix
end

@doc raw"""
    SizeDistribution

The type SizeDistribution abstracts the aerosol size distribution. The parameter A is 
a set of input parameters, e.g. for a lognormal function. The form contains a symbol 
that traces the function or process that created the distribution.

    A::Any                        # Input parameters [[N1,Dg1,σg1], ...] or DMA
    De::Vector{<:AbstractFloat}   # bin edges
    Dp::Vector{<:AbstractFloat}   # bin midpoints
    ΔlnD::Vector{<:AbstractFloat} # ΔlnD of the grid
    S::Vector{<:AbstractFloat}    # spectral density
    N::Vector{<:AbstractFloat}    # number concentration per bin
    form::Symbol                  # form of the size distribution [:lognormal, ....]

SizeDistributionscan be created by hand or through one of the constructor functions:
- [lognormal](@ref)
- [DMALognormalDistribution](@ref)
- [triangular](@ref)
"""
mutable struct SizeDistribution
    A::Any                        # Input parameters [[N1,Dg1,σg1], ...] or DMA
    De::Vector{<:AbstractFloat}   # bin edges
    Dp::Vector{<:AbstractFloat}   # bin midpoints
    ΔlnD::Vector{<:AbstractFloat} # ΔlnD of the grid
    S::Vector{<:AbstractFloat}    # spectral density
    N::Vector{<:AbstractFloat}    # number concentration per bin
    form::Symbol                  # form of the size distribution [:lognormal, ....]
end

@doc raw"""
    Regvars

The type Regvars abstracts the inversion setup, including the convolution matrix,
the idenity matrix, the response function (residual vector) and the initial guess.
The matrix 𝐀'𝐀 is stored as precomputed matrix to avoid recomputing it when evaluating
the derivatives in the l-curve search. Setting the number of BLAS threads is experimental.

    𝐀::Matrix{Float64}     # Convolution matrix
    𝐈::Matrix{Float64}      # Identity matrix
    B::Array{Float64}      # residual vector
    X₀::Array{Float64}     # initial guess
    AA::Matrix{Float64}    # precomputed 𝐀'𝐀 for speed
    n::Int                 # Blas threads

Regvars is initialized in the rinv function. See examples folder on how to use this 
structure at the top level.
"""
struct Regvars
    𝐀::Matrix{Float64}     # Convolution matrix
    𝐈::Matrix{Float64}      # Identity matrix
    B::Array{Float64}      # residual vector
    X₀::Array{Float64}     # initial guess
    AA::Matrix{Float64}    # precomputed A'A for speed
    n::Int                 # Blas # threads
end

### Constants
const kb = 1.38064852e-23    # Boltzmann constant      [J K-1 molec-1]
const r = 8.314472           # Universal gas const     [J K-1 mol-1]
const ec = 1.602176565e-19   # Elemental charge        [C]
const eps = 8.854187817e-12  # Dielectric constant     [farad m-1]
const na = 6.02214086e23     # Avagadro's constant     [molecule mol-1]
const mwair = 28.9647e-3     # Molecular weight of air [kg mol-1]
const di = 1e-7              # Initial guess diameter  [m]

# Coefficients for charge polynomials taken from TSI Manual.
# n2: negative 2 charges, p2: positive 2 charges on PARTICLE
const ax = DataFrame(
    n2 = [-26.3328, 35.9044, -21.4608, 7.0867, -1.3088, 0.1051],
    n1 = [-2.3197, 0.6175, 0.6201, -0.1105, -0.1260, 0.0297],
    p1 = [-2.3484, 0.6044, 0.4800, 0.0013, -0.1553, 0.0320],
    p2 = [-44.4756, 79.3772, -62.8900, 26.4492, -5.7480, 0.5049],
)

### Source files
include("dmafunctions.jl")
include("tdmafunctions.jl")
include("aerosolsizedistributions.jl")
include("regularization.jl")
include("loadtsidata.jl")
include("coagulationrates.jl")
include("benchmark.jl")
include("myfigure.jl")

end
