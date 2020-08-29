push!(LOAD_PATH, "../src/")
using Documenter, DifferentialMobilityAnalyzers, LinearAlgebra

makedocs(
  sitename = "DifferentialMobilityAnalyzers",
  authors = "Markus Petters",
  pages = Any[
    "Home" => "index.md",
    "Library" => Any[
      "Operators" => "operators.md",
      "Benchmarks" => "benchmarks.md"
    ]
  ]
)

