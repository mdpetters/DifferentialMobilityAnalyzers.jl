push!(LOAD_PATH, "../src/")
using Documenter, DifferentialMobilityAnalyzers, LinearAlgebra

makedocs(
  sitename = "DifferentialMobilityAnalyzers.jl",
  authors = "Markus Petters",
  pages = Any[
    "Home" => "index.md",
    "Getting Started" => Any[
      "Quick Start" => "quickstart.md",
      "Tutorial" => "tutorial.md",
      "Notebooks" => "notebooks.md"
    ],
    "Library" => Any[
      "Physics" => "physics.md",
      "Data Types" => "types.md",
      "Operators" => "operators.md",
      "Helper Functions" => "helpers.md",
      "Benchmarks" => "benchmarks.md"
    ]
  ]
)

