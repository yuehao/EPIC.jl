using Revise, EPIC
using Test
using BenchmarkTools
using StructArrays
using Plots
using LinearAlgebra

# test bunched beam with IP optics and 6D Gaussian distribution

pbeam=BunchedBeam(PROTON, 1e11, 250e9,  1000000, [2.5e-6, 1e-6, 1e-1])
pbeam.dist.x
a=ps6d{Float64}(1.0,2.0,3.0,4.0,5.0,6.0)
opIP=optics4DUC(4.0,0.0,1.0,0.0)
mainRF=AccelCavity(197e6, 1e6, 2520.0, 0.0)
initilize_6DGaussiandist!(pbeam, opIP, mainRF, 1.8e-3)
@benchmark begin
    initilize_6DGaussiandist!(pbeam, opIP, mainRF, 1.8e-3)
    #track!(pbeam.dist, mainRF, 1.8e-3)
end
get_2nd_moments!(pbeam)
get_centroid!(pbeam)
@benchmark begin get_emittance!(pbeam) end
pbeam.emittance
pbeam.moment2nd



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


## test transfer map
oneturn=TransferMap4D(opIP, 2*pi*0.08, 2*pi*0.139)
otm=oneturn.linearmap
otm[1,2]/otm[2,1]
@benchmark begin track!(pbeam, oneturn) end
@benchmark begin track!(pbeam.dist, oneturn) end

oneturnrad=OneTurnRadiation(4000.0, 4000.0, 2000.0)
track!(pbeam.dist, oneturnrad, [1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5])



# test Weak-Strong Simulation
@benchmark begin e,g=Bassetti_Erskine.(pbeam.dist.x,pbeam.dist.y,5.0,1.0) end
estrong=StrongThinGaussianBeam(1.0, 1.0, 0.9)

track!(pbeam.dist, estrong)

@testset "EPIC.jl" begin
    # Write your tests here.
    a=1
end
