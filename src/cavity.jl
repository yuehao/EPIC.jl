
import StaticArrays
struct CrabCavity <: AbstractCavity
    f::Float64
    v::Float64
    k::Float64
    ϕ::Float64
    errors::AbstractArray
    CrabCavity(freq,nv)=new(freq,nv,2*π*freq/2.99792458e8, phase, StaticArrays.@MVector [0.0,0.0])
end

function track(coor6d::AbstractVector, cc::CrabCavity)
    s,c=sincos(cc.k*coor6d[5]+cc.ϕ)
    coor6d[2]+= cc.v*s    #px
    coor6d[6]+= cc.v*coor6d[1]*c  # delta
end


