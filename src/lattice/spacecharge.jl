using SpecialFunctions
using Statistics
using StaticArrays

struct SpaceChargeLens <: AbstractSpaceCharge
    optics::AbstractOptics4D   # optics at the location of the lens
    ds::Float64  # effective length of the lens
    SpaceChargeLens(optics, ds) = new(optics, ds)
end

function track!(dist::AbstractVector{ps6d{T}}, σx::Float64, σy::Float64, factor::Float64, σz::Float64, temp1) where T
    
    fieldvec_thread=[MVector{3}(0.0, 0.0, 0.0)  for j = 1:Threads.nthreads()]


    @inbounds Threads.@threads :static for j in eachindex(dist.x)

        temp1[j] = 1.0 /(σz * sqrt(2*pi)) * exp((-0.5) * dist.z[j]^2 /σz^2) # lambda_z
            

        Bassetti_Erskine!(fieldvec_thread[Threads.threadid()], dist.x[j], dist.y[j], σx , σy) #function defined in strongbb.jl

        # gaussian function for part. distribution, lambda_z  (bbSC.dist.z[j]= temp1[j])
        
        
        #factor1 = 2*bbSC.particle.classrad0/(bbSC.beta^2* bbSC.gamma^3)
        

        # delta p_/p_0
        dist.px[j] += factor* temp1[j]* fieldvec_thread[Threads.threadid()][1]
        dist.py[j] += factor* temp1[j]* fieldvec_thread[Threads.threadid()][2]


    end



    
    return nothing
end

function track!(beam::BunchedBeam, sc::SpaceChargeLens)
    get_emittance!(beam)
    σx = sqrt(beam.emittance[1]*sc.optics.optics_x.beta)
    σy = sqrt(beam.emittance[2]*sc.optics.optics_y.beta)
    σz = beam.beamsize[5]

    factor=beam.particle.classrad0/beam.gamma^3/beam.beta^2*sc.ds*beam.num_particle
    track!(beam.dist, σx, σy, factor, σz, beam.temp1)
end