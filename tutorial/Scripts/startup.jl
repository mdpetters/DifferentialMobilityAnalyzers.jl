# This startup file loads the required packages for the tutorial and defines
# the display function to embed videos into Jupyter using Julia

using WebIO
using PyCall
using DifferentialMobilityAnalyzers
using Printf
using LsqFit
using Gadfly
using LinearAlgebra
using DataFrames
using Random
using Distributions
using SpecialFunctions
using JLD2
using FileIO
using CSV
using Interpolations
using MLStyle

WebIO.setup(:ijulia)
IPyDisplay = pyimport("IPython.display")
play(audiofile::AbstractString) = IPyDisplay.Audio(audiofile)

display(url::AbstractString; res = :240) = @match res begin
    :240 => IPyDisplay.YouTubeVideo((split(url, "="))[2], width = 426, height = 240)
    :360 => IPyDisplay.YouTubeVideo((split(url, "="))[2], width = 640, height = 360)
    :480 => IPyDisplay.YouTubeVideo((split(url, "="))[2], width = 824, height = 480)
end

matrices = JLD2.jldopen("precomputed_matrices.jld", "r") do file
    read(file, "matrices")
end
