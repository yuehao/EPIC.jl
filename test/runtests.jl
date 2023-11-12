using EPIC
using Test
using BenchmarkTools
using Revise
using Plots

# test bunched beam with IP optics and 6D Gaussian distribution

pbeam=BunchedBeam(PROTON, 1e11, 250e9,  400000, [4e-6, 1e-6, 1e-1])
opIP=optics4DUC(4.0,0.0,1.0,0.0)
mainRF=AccelCavity(197e6, 1e6, 2520.0, 0.0)
@btime begin mat=initilize_6DGaussiandist!(pbeam, opIP, mainRF, 1.8e-3) end
sqrt(sum(pbeam.dist.dp .* pbeam.dist.dp)/pbeam.num_macro)
sigz=sqrt(sum(pbeam.dist.z .* pbeam.dist.z)/pbeam.num_macro)

# test Strong Beam setup
pstbeam=StrongGaussianBeam(PROTON, 1e11, 250e9,  opIP, [4e-6, 1e-6, sigz], 20)
initilize_zslice!(pstbeam, :gaussian, :evennpar, 100.0)
pstbeam.zslice_center


initilize_zslice!(pstbeam, pbeam.dist.z, :evennpar)
pstbeam.zslice_npar
pstbeam.zslice_center

# test histogram and wakefield
@btime begin
    zhist, zhist_edges=histogram1DinZ!(pbeam, 200)
end
plot(zhist_edges[1:end-1], zhist, seriestype=:scatter, markersize=1, legend=false, xlabel="z [m]", ylabel="N", title="6D Gaussian distribution")

RLCwake = LongitudinalRLCWake(1e9, 1e5, 1.0)
old_dp=pbeam.dist.dp.*1.0
track!(pbeam, RLCwake)

plot(pbeam.dist.z, pbeam.dist.dp.-old_dp, markersize=1, seriestype=:scatter, legend=false, xlabel="z [m]", ylabel="dp/p" )


@testset "EPIC.jl" begin
    # Write your tests here.
    a=1
end
