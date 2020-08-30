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
      "Size Distributions" => "man/sizedistribution.md",
      "Initializing DMAs" => "man/dmas.md",
      "Tranmsission through the DMA" => "man/transmission.md",
      "Creating Convolution Matrices" => "man/matrix.md",
      "Size Distribution Inversion" => "man/inversion.md",
      "Tandem DMA Setups" => "man/tandem.md"
    ],
    "Library" => Any[
      "Physical Relations" => "lib/physics.md",
      "Data Types" => "lib/types.md",
      "Operators" => "lib/operators.md",
      "Helper Functions" => "lib/helpers.md",
      "Inversion Routines" => "lib/inversion.md",
      "Benchmarks" => "lib/benchmarks.md"
    ]
  ]
)

