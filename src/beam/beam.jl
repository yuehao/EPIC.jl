
include("particle.jl")
include("optics.jl")
include("dist.jl")

abstract type AbstractBeam end

struct StrongBeam <: AbstractBeam
    particle::ParticleType
    num_particle::Float64
    total_energy::Float64
    momentum::Float64
    gamma::Float64
    beta::Float64
    function StrongBeam(particle::ParticleType, np::Float64, energy::Float64)
        momentum=sqrt(energy*energy-particle.mass*particle.mass)
        gamma=energy/particle.mass
        beta=momentum/energy
        new(particle,np,energy,momentum,gamma,beta)
    end
end

mutable struct  BunchedBeam <: AbstractBeam
    particle::ParticleType
    num_particle::Float64
    total_energy::Float64
    momentum::Float64
    gamma::Float64
    beta::Float64
    num_macro::Int64
    dist::Array{Float64, 2}
    function BunchedBeam(particle::ParticleType, np::Float64, energy::Float64, nmacro::Int)
        momentum=sqrt(energy*energy-particle.mass*particle.mass)
        gamma=energy/particle.mass
        beta=momentum/energy
        new(particle,np, energy,momentum,gamma,beta,nmacro,zeros(6,nmacro))
    end
end



function initilize_6ddist!(6ddist::Array{Float64, 2}, opxy::Optics2DUncoupled, opz::Optics1D)
end

function initilize_6ddist!(6ddist::Array{Float64, 2}, opxy::Optics2DUncoupled, rf::AccelCavity)
end