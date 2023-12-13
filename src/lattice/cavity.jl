

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
    ϕs::Float64 # synchronous phase π/2 for accelerating on crest
    
    AccelCavity(freq,nv,h,ϕs)=new(freq,nv, 2*π*freq/2.99792458e8, h, ϕs)
end


function track!(ps6dcoor::AbstractVector, cc::CrabCavity)
   sv,cv .= sincos.(cc.k .* ps6dcoor.z .+ cc.ϕ)
   ps6dcoor.px .+= cc.v .* sv    #px
   ps6dcoor.dp += cc.v .* ps6dcoor.x .* cv  # delta
end

function track!(ps6dcoor::AbstractVector{ps6d{T}}, ac::AccelCavity, b2E::Float64) where T
   v_β2E=ac.v/b2E
   sv = sin.(ac.k .* ps6dcoor.z .+ ac.ϕs) .- sin(ac.ϕs)
   ps6dcoor.dp .+= v_β2E .* sv   # update delta p
end

function track!(beam::AbstractBeam, cc::CrabCavity)
   s,c=sincos.(cc.k*beam.dist.z+cc.ϕ)
   beam.dist.px .+= cc.v .* s    #px
   beam.dist.dp .+= cc.v .* beam.dist.x .* c  # delta
end

function track!(beam::AbstractBeam, ac::AccelCavity)
   v_β2E=ac.v/(beam.beta*beam.beta*beam.total_energy)
   ss=sin.(ac.k .* beam.dist.z .+ ac.ϕs) .- sin(ac.ϕs)
   beam.dist.dp .+= v_β2E .* ss   # update delta p
end



