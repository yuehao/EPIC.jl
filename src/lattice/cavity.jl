using StaticArrays

mutable struct CrabCavity <: AbstractCavity
    f::Float64
    v::Float64
    k::Float64
    ϕ::Float64
    errors::AbstractArray # 1: Voltage error, 2: Phase error
    CrabCavity(freq,nv)=new(freq,nv,2*π*freq/2.99792458e8, 0.0, StaticArrays.@MVector [0.0,0.0])
end

mutable struct AccelCavity <: AbstractCavity
    f::Float64
    v::Float64
    k::Float64
    h::Int64
    ϕs::Float64
    
    AccelCavity(freq,nv,h,ϕs)=new(freq,nv, 2*π*freq/2.99792458e8, h, ϕs)
end



function track(coor6d::AbstractVector, cc::CrabCavity)
    s,c=sincos(cc.k*coor6d[5]+cc.ϕ)
    coor6d[2]+= cc.v*s    #px
    coor6d[6]+= cc.v*coor6d[1]*c  # delta
end

function track(coor6d::AbstractVector, ac::AccelCavity, β2E::Float64)
    ss=sin(ac.k*coor6d[5]+ac.ϕs)-sin(ac.ϕs)
    coor6d[6]+= ac.v*ss / β2E   #px
end



