using BenchmarkTools

"""
    benchmark(bins::Integer, num_threads::Integer)

Computes benchmarks for the three slowest operations, rinv, setupDMA, setupSMPS
- bins is the number of DMA bins
- num_threads is the number of BLAS threads

The function returns a dataframe that includes cpuinfo, juliaversion, blasvendor,
blasthreads, number of bins, and the three timed benchmarks for rinv, setupDMA, setupSMPS.
"""
function benchmark(bins::Integer, num_threads::Integer)
    # Load a simple comma delimited text file
    path = @__DIR__
    df = CSV.read(path*"/example_data.csv", DataFrame)

    # Setup the DMA
    t, p, lpm = 293.15, 940e2, 1.666e-5      
    r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369   
    Λ = DMAconfig(t, p, 1lpm, 4lpm, r₁, r₂, l, 0.0, :+, 6, :cylindrical)
    δ = setupDMA(Λ, vtoz(Λ, 10000), vtoz(Λ, 10), bins)
    
    𝕣 = (df, :Dp, :Rcn, δ) |> interpolateDataFrameOntoδ

    a = @benchmark setupDMA($Λ, vtoz($Λ, 10000), vtoz($Λ, 10), $bins)
    b = @benchmark setupSMPS($Λ, 10000, 10, $bins, 1.0)
    c = @benchmark rinv($(𝕣.N), $δ; λ₁ = 0.1, λ₂ = 1.0)
    d = @benchmark rinv2($(𝕣.N), $δ; λ₁ = 0.1, λ₂ = 1.0)
 
    cpuio = IOBuffer() # print cpu_summary with correct alignment
    Sys.cpu_summary(cpuio)
    CPU = String((split(String(take!(cpuio)), "\n"))[1])

    juliav = "v$VERSION"
    blasvendor = String(BLAS.vendor())

    DataFrame(
        cpuinfo = CPU,
        julia = juliav,
        blas = blasvendor,
        threads = num_threads,
        nbins = bins,
        # setupDMA = string(BenchmarkTools.prettytime(median(a).time)),
        # setupSMPS = string(BenchmarkTools.prettytime(median(b).time)),
        rinv = string(BenchmarkTools.prettytime(median(c).time)),
        rinv2 = string(BenchmarkTools.prettytime(median(d).time))
    )
end

"""
    runbenchmarks()

Runs a set of standard benchmarks, varying the number bins and returns a dataframe with
results. This function may take several minutes to complete. 
"""
function runbenchmarks()
    try 
        benchmark(10,1)
    catch
    end
    bins = [30, 60, 120]
    xx = map(n->benchmark(n,1), bins)
    return vcat(xx...)
end
