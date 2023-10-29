abstract type AbstractStrongBeamBeam <:AbstractElement end
using SpecialFunctions
using Statistics

struct StrongThinGaussianBeam <: AbstractStrongBeamBeam
    amplitude::Float64
    rmssizex::Float64
    rmssizey::Float64
    zloc::Float64
    xoffset::Float64
    yoffset::Float64
    StrongThinGaussianBeam(amp::Float64, rx::Float64, ry::Float64, zloc::Float64=0.0, xoff::Float64=0.0, yoff::Float64=0.0)=new(amp,rx,ry,zloc,xoff,yoff)
end

struct StrongGaussianBeam <: AbstractStrongBeamBeam  # Strong Beam with transverse Gaussian distribution
    particle::ParticleType
    num_particle::Float64  # Number of particles
    total_energy::Float64 # Total energy of the beam
    momentum::Float64  # Design Momentum of the beam
    gamma::Float64  # Relativistic gamma
    beta::Float64  # Relativistic beta v/c
    optics::AbstractOptics4D # optics @IP
    beamsize::Vector{Float64} # Beam size at IP
    nzslice::Int64 # Number of slices in z direction
    zslice_center::Vector{Float64} # z center of each slice
    zslice_npar::Vector{Float64} # amplitude of each slice
    function StrongGaussianBeam(particle::ParticleType, np::Float64, energy::Float64, op::AbstractOptics4D, bs::Vector{Float64}, nz::Int)
        momentum=sqrt(energy*energy-particle.mass*particle.mass)
        gamma=energy/particle.mass
        beta=momentum/energy
        new(particle,np,energy,momentum,gamma,beta, op, bs, Int64(nz), zeros(nz), zeros(nz))
    end  
end

function initilize_zslice!(beam::StrongGaussianBeam, profile::Symbol, slice_type::Symbol, zrange::Float64=5.0)
    zmin=-zrange*beam.beamsize[3]
    zmax=zrange*beam.beamsize[3]
    if profile == :gaussian
        if slice_type == :evenzsep
            zedge=collect(range(zmin, stop=zmax, length=beam.nzslice+1))
            beam.zslice_center .= 0.5.*(zedge[1:end-1]+zedge[2:end])
            beam.zslice_npar .= exp.(-0.5.*(beam.zslice_center.^2)./beam.beamsize[3]^2)
            beam.zslice_npar .= beam.zslice_npar./sum(beam.zslice_npar).*beam.num_particle
        end
        if slice_type == :evennpar
            npartedge=collect(range(0.0, stop=1.0, length=beam.nzslice+1))
            npartcenter=0.5.*(npartedge[1:end-1]+npartedge[2:end])  
            beam.zslice_center .= (sqrt(2.0)*beam.beamsize[3]).*erfinv.(2.0.*npartcenter.-1.0)
            beam.zslice_npar .= zeros(beam.nzslice).+1.0/beam.nzslice*beam.num_particle
        end
    end

    if profile == :uniform
        zedge=collect(range(zmin, stop=zmax, length=beam.nzslice+1))
        beam.zslice_center .= 0.5.*(zedge[1:end-1]+zedge[2:end])
        beam.zslice_npar .= zeros(beam.nzslice).+1.0/beam.nzslice*beam.num_particle
    end
    
end

function initilize_zslice!(beam::StrongGaussianBeam, zlist::Vector{Float64}, slice_type::Symbol)
    zmin=minimum(zlist)
    zmax=maximum(zlist)
    if slice_type == :evenzsep
        zedge=collect(range(zmin, stop=zmax, length=beam.nzslice+1))
        beam.zslice_center=0.5.*(zedge[1:end-1]+zedge[2:end])
        beam.zslice_npar=zeros(beam.nzslice)
        for i in 1:beam.nzslice
            beam.zslice_npar[i]=sum((zlist.>=zedge[i]).&(zlist.<zedge[i+1]))
        end
        beam.zslice_npar=beam.zslice_npar./sum(beam.zslice_npar).*beam.num_particle
    end
    if slice_type == :evennpar
        sort_zlist=sort(zlist)
        nzlist=length(sort_zlist)
        npartedge=Int.(floor.(collect(range(0, nzlist, length=beam.nzslice+1))))
        beam.zslice_npar .= (npartedge[2:end]-npartedge[1:end-1])/nzlist*beam.num_particle
        beam.zslice_center .= zeros(beam.nzslice)
        for i in 1:beam.nzslice
            beam.zslice_center[i]=mean(sort_zlist[npartedge[i]+1:npartedge[i+1]])
        end
    end
end


function Bassetti_Erskine(x::Float64, y::Float64, σx::Float64, σy::Float64)
	sqrtδσ2=sqrt(Complex(2*(σx*σx-σy*σy)))
	term1=erfcx(-1im*(x+1im*y)/sqrtδσ2)
	term2=erfcx(-1im*(x*σy/σx+1im*y*σx/σy)/sqrtδσ2)
	termexp=exp(-x*x/2/σx/σx-y*y/2/σy/σy)
	complex_e=-1im*2*sqrt(pi)/sqrtδσ2*(term1-termexp*term2)
	return real(complex_e), -imag(complex_e)
end

function track(coor6d::AbstractArray, stgb::StrongThinGaussianBeam)
    sloc=(stgb.zloc .+ coor6d.z)/2.0
    coor6d.x .+= (coor6d.px .* stgb.zloc)
    coor6d.y .+= (coor6d.py .* stgb.zloc)
    ex,ey=Bassetti_Erskine.(coor6d.x.-stgb.xoffset, coor6d.y.-stgb.yoffset, stgb.rmssizex, stgb.rmssizey)
    coor6d.px .+= stgb.amplitude .* ex
    coor6d.py .+= stgb.amplitude .* ey
    coor6d.x .-= (coor6d.px .* stgb.zloc)
    coor6d.y .-= (coor6d.py .* stgb.zloc)
end




