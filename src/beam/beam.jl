


include("particle.jl")


struct BunchedBeam<: AbstractBeam 
    particle::ParticleType
    num_particle::Float64
    total_energy::Float64  
    momentum::Float64  # Design Momentum of the beam
    gamma::Float64  # Relativistic gamma
    beta::Float64  # Relativistic beta v/c
    num_macro::Int64 # Number of macro particles
    dist::StructArray{ps6d{Float64}} # 6D distribution
    inzindex::Vector{Int64} # in which z slice

    emittance::Vector{Float64} # emittance in x, y, z
    centroid::Vector{Float64} # centroid in x, px, y, py, z, pz
    moment2nd::Matrix{Float64} # 2nd moment matrix

    function BunchedBeam(particle::ParticleType, np::Float64, energy::Float64, nmacro::Int, 
                        emittance::Vector{Float64}, centroid::Vector{Float64}=zeros(6))
        momentum=sqrt(energy*energy-particle.mass*particle.mass)
        gamma=energy/particle.mass
        beta=momentum/energy
        new(particle, np, energy,momentum,gamma,beta,nmacro,
            StructArray{ps6d{Float64}}(undef,nmacro), 
            zeros(Int64,nmacro), emittance, centroid, zeros(6,6))
    end
end


include("dist.jl")


