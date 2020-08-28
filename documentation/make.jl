push!(LOAD_PATH, "../src/")
using Documenter


makedocs(
  sitename = "DifferentialMobilityAnalyzers.jl",
  authors = "Markus Petters",
  pages = [
    "Home" => "index.md",
  ]

)

