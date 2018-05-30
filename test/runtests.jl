using DifferentialMobilityAnalyzers, Distributions
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

tests = ["dmafunctions",
         "sizedistribution",
         "inversion",
         "coagulation",
         "miscellaneous"]

println("Running tests:")

for t in tests
    println(" * $(t)")
    include("$(t).jl")
end
