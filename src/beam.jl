import .optics

abstract type AbstractBeam end

struct  BunchBeam <: AbstractBeam
    num_particle::Float64
    mass::Float64
    total_energy::Float64
    momentum::Float64
    gamma::Float64
    beta::Float64
    
    num_macro::Int64
    dist::Array{Float64, 2}
    function BunchBeam(np, mass, energy, nm)
        momentum=sqrt(energy*energy-mass*mass)
        gamma=energy/mass
        beta=momentum/energy
        new(np,mass,energy,momentum,gamma,beta,nm,zeros(6,nm))
    end
end



