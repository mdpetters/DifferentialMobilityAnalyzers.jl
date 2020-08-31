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
    "Manual" => Any[
      "Overview" => "man/overview.md",
      "Size Distributions" => "man/sizedistribution.md",
      "Initializing DMAs" => "man/dmas.md",
      "Convolution Matrices" => "man/matrix.md",
      "Forward Models" => "man/forward.md",
      "Inverse Models" => "man/inverse.md",
      "Performance Benchmarks" => "man/benchmarks.md",
    ],
    "Library" => Any[
      "Physics" => "lib/physics.md",
      "Data Types" => "lib/types.md",
      "Operators" => "lib/operators.md",
      "Helper Functions" => "lib/helpers.md",
      "Inversion Routines" => "lib/inversion.md",
      "Benchmarks" => "lib/benchmarks.md"
    ],
    "References" => "references.md"
  ]
)

