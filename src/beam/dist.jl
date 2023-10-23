using LinearAlgebra


const distribution1D_type=Set([:custom, :gaussian, :uniform])
const zslice_type=Set([:evenzsep, :evennpar])


function initilize_zslice!(beam::StrongGaussianBeam, profile::Symbol, slice_type::Symbol, zrange::Float64=5.0)
    zmin=-zrange*beam.beamsize[3]
    zmax=zrange*beam.beamsize[3]
    if profile == :gaussian
        if slice_type == :evenzsep
            zedge=collect(range(zmin, stop=zmax, length=beam.nzslice+1))
            beam.zslice_center=0.5.*(zedge[1:end-1]+zedge[2:end])
            beam.zslice_npar=exp.(-0.5.*(beam.zslicecenter.^2)./beam.beamsize[3]^2)
            beam.zslice_npar=beam.zslice./sum(beam.zslice).*beam.num_particle
        end
        if slice_type == :evennpar
            npartedge=collect(range(0.0, stop=1.0, length=beam.nzslice+1))
            npartcenter=0.5.*(npartedge[1:end-1]+npartedge[2:end])  
            beam.zslicecenter=(sqrt(2.0)*beam.beamsize[3]).*erfinv.(2.0.*npartcenter.-1.0)
            beam.zslice_npar=zeros(beam.nzlice).+1.0/beam.nzslice*beam.num_particle
        end
    end

    if profile == :uniform
        zedge=collect(range(zmin, stop=zmax, length=beam.nzslice+1))
        beam.zslice_center=0.5.*(zedge[1:end-1]+zedge[2:end])
        beam.zslice_npar=zeros(beam.nzslice).+1.0/beam.nzslice*beam.num_particle
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
        beam.zslice_npar=(npartedge[2:end]-npartedge[1:end-1])/nzlist*beam.num_particle
        beam.zslice_center=zeros(beam.nzslice)
        for i in 1:beam.nzslice
            beam.zslice_center[i]=mean(sort_zlist[npartedge[i]+1:npartedge[i+1]])
        end
    end
end

strongebeam=StrongGaussianBeam(ELECTRON, 1e10, 6500.0e9, optics4DUC(0.2, 0.05, 0.1, 0.), [1.0, 1.0, 1.0], 100)

function initilize_6DGaussiandist!(beam::BunchedBeam, optics::optics4DUC, rf::AccelCavity, αc::Float64)
    # 6D Gaussian distribution
    beam.dist[1,:]=randn(beam.num_macro)
    beam.dist[2,:]=randn(beam.num_macro)
    beam.dist[3,:]=randn(beam.num_macro)
    beam.dist[4,:]=randn(beam.num_macro)
    beam.dist[5,:]=randn(beam.num_macro)
    beam.dist[6,:]=randn(beam.num_macro)

    moment_matrix=zeros(6,6)
    for i in 1:6
        for j in i:6
            moment_matrix[i,j]=sum(beam.dist[i,:].*beam.dist[j,:])/beam.num_macro
        end
    end
    moment_matrix=moment_matrix+moment_matrix'-diagm(diag(moment_matrix))
    eigval,eigvec=eigen(moment_matrix)
    mscale=eigvec * diagm(sqrt.(eigval)) * eigvec'
    beam.dist=mscale*beam.dist

    beam.dist[1,:]=beam.dist[1,:].*sqrt(beam.emittance[1]*optics.optics_x.beta)
    beam.dist[2,:]=beam.dist[2,:].*sqrt(beam.emittance[1]/optics.optics_x.beta)
    beam.dist[2,:] += beam.dist[1,:].*(optics.optics_x.alpha/optics.optics_x.beta)

    beam.dist[3,:]=beam.dist[3,:].*sqrt(beam.emittance[2]*optics.optics_y.beta)
    beam.dist[4,:]=beam.dist[4,:].*sqrt(beam.emittance[2]/optics.optics_y.beta)
    beam.dist[4,:] += beam.dist[3,:].*(optics.optics_y.alpha/optics.optics_y.beta)

    #generate longtiudinal distribution based on small amplitude approximation
    eta_p=αc-1.0/beam.gamma^2
    ωrf=2.0*π*rf.f
    Qs=sqrt(rf.v*rf.h*abs(eta_p*cos(rf.ϕs))/2/π/beam.beta0^2/beam.energy)

    emit_deltap_z=beam.emittance[3]*2.99792458e8/beam.beta0/beam.energy
    invbeta_deltap_z=Qs*rf.k/rf.h/abs(eta_p)
    beam.dist[5,:]=beam.dist[5,:].*sqrt(emit_deltap_z/invbeta_deltap_z)
    beam.dist[6,:]=beam.dist[6,:].*sqrt(emit_deltap_z*invbeta_deltap_z)
    
    
end
