abstract type AbstractCavity <:AbstractElement end

mutable struct CrabCavity <: AbstractCavity
    f::Float64 # frequency
    v::Float64 # voltage
    k::Float64 # wavenumber
    ϕ::Float64 # phase
    errors::AbstractArray # 1: Voltage error, 2: Phase error
    CrabCavity(freq,nv)=new(freq,nv,2*π*freq/2.99792458e8, 0.0, StaticArrays.@MVector [0.0,0.0])
end

mutable struct AccelCavity <: AbstractCavity
    f::Float64 # frequency
    v::Float64 # voltage
    k::Float64 # wavenumber
    h::Int64 # harmonic number
    ϕs::Float64 # synchronous phase
    
    AccelCavity(freq,nv,h,ϕs)=new(freq,nv, 2*π*freq/2.99792458e8, h, ϕs)
end



function track(coor6d::AbstractArray, cc::CrabCavity)
    s,c=sincos(cc.k*coor6d.z+cc.ϕ)
    coor6d.px+= cc.v*s    #px
    coor6d.dp+= cc.v*coor6d[1]*c  # delta
end

function track(coor6d::AbstractArray, ac::AccelCavity, β2E::Float64)
    ss=sin(ac.k*coor6d.z+ac.ϕs)-sin(ac.ϕs)
    coor6d.dp+= ac.v*ss / β2E   # update delta p
end



