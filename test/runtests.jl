using DifferentialMobilityAnalyzers, Distributions, Test

tests = [
    "dmafunctions",
    "sizedistribution",
    "inversion1",
    "inversion2",
    "inversion3",
    "coagulation",
]

println("Running tests:")

for testfun in tests
    println(" * $(testfun)")
    include("$(testfun).jl")
end
