using EPIC
using Test
using StructArrays
using BenchmarkTools

# test bunched beam with IP optics and 6D Gaussian distribution

pbeam=BunchedBeam(PROTON, 1e11, 250e9,  1000000, [4e-6, 1e-6, 1e-1])
pbeam.emittance

opIP=optics4DUC(4.0,0.0,1.0,0.0)
mainRF=AccelCavity(197e6, 1e6, 2520.0, 0.0)
mat=initilize_6DGaussiandist!(pbeam, opIP, mainRF, 1.8e-3)

sqrt(sum(pbeam.dist.dp .* pbeam.dist.dp)/pbeam.num_macro)
sigz=sqrt(sum(pbeam.dist.z .* pbeam.dist.z)/pbeam.num_macro)

# test Strong Beam setup
pstbeam=StrongGaussianBeam(PROTON, 1e11, 250e9,  opIP, [4e-6, 1e-6, sigz], 20)
initilize_zslice!(pstbeam, :gaussian, :evennpar, 100.0)
pstbeam.zslice_center


initilize_zslice!(pstbeam, pbeam.dist.z, :evennpar)
pstbeam.zslice_npar
pstbeam.zslice_center


using StructArrays
using BenchmarkTools
using StaticArrays
n=1000000
a=StructArray{ps6d}((randn(n),randn(n),randn(n),randn(n),randn(n),randn(n)))  
dim=:x
k=:(a.$dim)
eval(k)
@btime begin
    b=Matrix{Float64}(undef,6,n)
    b[1,:]=a.x
    b[2,:]=a.px
    b[3,:]=a.y
    b[4,:]=a.py
    b[5,:]=a.z
    b[6,:]=a.dp
end


@btime begin
    b=Matrix{Float64}(undef,n,6)
    b[:,1]=a.x
    b[:,2]=a.px
    b[:,3]=a.y
    b[:,4]=a.py
    b[:,5]=a.z
    b[:,6]=a.dp
end


@testset "EPIC.jl" begin
    # Write your tests here.
    a=1
end
