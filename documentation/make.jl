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
      "Initializing DMAs" => "dmas.md",
      "Tranmsission through the DMA" => "transmission.md",
      "Size Distribution Inversion" => "inversion.md",
    ],
    "Library" => Any[
      "Physical Relations" => "physics.md",
      "Data Types" => "types.md",
      "Operators" => "operators.md",
      "Helper Functions" => "helpers.md",
      "Benchmarks" => "benchmarks.md"
    ]
  ]
)

