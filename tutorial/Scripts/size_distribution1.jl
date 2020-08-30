r = x -> Int.(round.(x, digits = 0))  # Function to round and convert to Int
𝕟 = lognormal([[200, 80, 1.2]]; d1 = 30.0, d2 = 300.0, bins = 10);
DataFrame(
    Dlow = r(𝕟.De[1:end-1]),
    Dup = r(𝕟.De[2:end]),
    ΔlnD = round.(𝕟.ΔlnD, digits = 2),
    Dp = r(𝕟.Dp),
    S = r(𝕟.S),
    N = r(𝕟.N),
)
