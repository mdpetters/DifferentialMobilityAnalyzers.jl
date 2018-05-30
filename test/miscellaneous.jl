pwd()
cd("../docs/")
𝕟, rawdp, rawc, Nt_TSI = loadtsidata();
@test round.(Int,Nt_TSI[[1,10]]) == [197,4979]
cd("../test/")
