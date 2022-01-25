using DifferentialMobilityAnalyzers

function initializeDMAs(Dd, k)
	t, p = 295.15, 1e5
	qsa, qsh = 1.66e-5, 8.33e-5
	r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369
	Λ₁ = DMAconfig(t, p, qsa, qsh, r₁, r₂, l, 160.0, :-, 6, :cylindrical)
    Λ₂ = DMAconfig(t, p, qsa, qsh, r₁, r₂, l, 0.0, :-, 6, :cylindrical)
    δ₁ = setupDMA(Λ₁, dtoz(Λ₁, 500e-9), dtoz(Λ₁, 8e-9), 120)   
	δ₂ = setupDMA(Λ₂, dtoz(Λ₂, 5.0*Dd), dtoz(Λ₂, 0.8*Dd), k)
    Λ₁, Λ₂, δ₁, δ₂ 
end

function TDMAmatrix(𝕟ᶜⁿ, Dd, Λ₁, Λ₂, δ₂, bins)
	ge = δ₂.De./(Dd*1e9)
	Δgf = ge[1:end-1].-ge[2:end];
	gf = δ₂.Dp./(Dd*1e9)

	𝐀 = @>> begin
		TDMA1Ddomainfunction(𝕟ᶜⁿ, Λ₁, Λ₂, (Dd, 0.8, 5.0, bins));
	 	designmatrix(gf)
    end
    gf, ge, 𝐀
end

function poisson_noise(Qcpc, N; seed = seed, t = 2.0)
    Random.seed!(seed)
	lpm = 16.666666   # 16 cm3 s-1 = 1 lpm
	Q = Qcpc * lpm    # Flow in cm3 s-1
    c = N * Q * t   # number of counts in each bin
	map(c) do i
		f = rand(Poisson(i), 1)
		f[1] / (Q * t)
	end
end

Normalize(x) = x./sum(x)
RMSE(x,y) = sqrt(sum((x .- y).^2.0)./length(x))

function test_cases(case::String, gf, dg, bins)
	f = @match case begin
		"Bimodal"          => 0.7*pdf(Normal(1.3,0.07), gf) + pdf(Normal(1.7,0.2), gf)
		"Truncated Normal" => pdf(truncated(Normal(1.2,0.2) , 1, 17), gf)
		"Uniform"          => begin
			j = argmin(abs.(gf .- 1.3))
			i = argmin(abs.(gf .- 1.8))
			a = @> zeros(bins) setindex!(ones(j-i+1), i:j)
			a .* dg
		end
		"Single Channel"   => @> zeros(bins) setindex!(1.0, argmin(abs.(gf .- 1.6)))	
		"Two Channel"      => begin
			i = argmin(abs.(gf .- 1.6))
			a = @> zeros(bins) setindex!([0.7, 0.3], [i,i - 3])
			a .* dg
		end
	end
	Normalize(f)
end

function df_to_stack(mdf, Dlow, Dup)
    df = mdf[4:end, :]
    Dp = mdf[2, 5:end] |> Vector
    ii = [BitArray([0, 0, 0, 0]); (Dp .>= Dlow) .& (Dp .< Dup)]
    Dl = mdf[1, ii] |> Vector
    Dp = mdf[2, ii] |> Vector
    Du = mdf[3, ii] |> Vector

    n = length(Dp)
    mapreduce(vcat, 1:length(df[!, 1])-1) do i
        N = df[i, ii] |> Vector
        DataFrame(
            x1 = df[i, :timeInt64],
            x2 = df[i+1, :timeInt64],
            y1 = [Du[j] for j = 1:n],
            y2 = [Dl[j] for j = 1:n],
            N = maximum(N) > 0 ? N ./ maximum(N) : 0,
        )
    end
end
