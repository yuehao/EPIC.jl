using LinearAlgebra


const distribution1D_type=Set([:custom, :gaussian, :uniform])
const zslice_type=Set([:evenzsep, :evennpar])

function get_centroid!(beam::BunchedBeam)
    # calculate centroid
    beam.centroid[1] = sum(beam.dist.x)/beam.num_macro
    beam.centroid[2] = sum(beam.dist.px)/beam.num_macro
    beam.centroid[3] = sum(beam.dist.y)/beam.num_macro
    beam.centroid[4] = sum(beam.dist.py)/beam.num_macro
    beam.centroid[5] = sum(beam.dist.z)/beam.num_macro
    beam.centroid[6] = sum(beam.dist.dp)/beam.num_macro
end

function get_2nd_moments!(beam::BunchedBeam)
    # calculate emittance
    
    ## Better way to do this?  Really messy but works...
    beam.moment2nd .= 0.0
    beam.moment2nd[1,1] = sum(beam.dist.x .* beam.dist.x)/beam.num_macro 
    beam.moment2nd[1,2] = sum(beam.dist.x .* beam.dist.px)/beam.num_macro 
    beam.moment2nd[1,3] = sum(beam.dist.x .* beam.dist.y)/beam.num_macro 
    beam.moment2nd[1,4] = sum(beam.dist.x .* beam.dist.py)/beam.num_macro 
    beam.moment2nd[1,5] = sum(beam.dist.x .* beam.dist.z)/beam.num_macro 
    beam.moment2nd[1,6] = sum(beam.dist.x .* beam.dist.dp)/beam.num_macro 

    beam.moment2nd[2,2] = sum(beam.dist.px .* beam.dist.px)/beam.num_macro 
    beam.moment2nd[2,3] = sum(beam.dist.px .* beam.dist.y)/beam.num_macro 
    beam.moment2nd[2,4] = sum(beam.dist.px .* beam.dist.py)/beam.num_macro 
    beam.moment2nd[2,5] = sum(beam.dist.px .* beam.dist.z)/beam.num_macro 
    beam.moment2nd[2,6] = sum(beam.dist.px .* beam.dist.dp)/beam.num_macro 

    beam.moment2nd[3,3] = sum(beam.dist.y .* beam.dist.y)/beam.num_macro 
    beam.moment2nd[3,4] = sum(beam.dist.y .* beam.dist.py)/beam.num_macro 
    beam.moment2nd[3,5] = sum(beam.dist.y .* beam.dist.z)/beam.num_macro 
    beam.moment2nd[3,6] = sum(beam.dist.y .* beam.dist.dp)/beam.num_macro 

    beam.moment2nd[4,4] = sum(beam.dist.py .* beam.dist.py)/beam.num_macro 
    beam.moment2nd[4,5] = sum(beam.dist.py .* beam.dist.z)/beam.num_macro 
    beam.moment2nd[4,6] = sum(beam.dist.py .* beam.dist.dp)/beam.num_macro 

    beam.moment2nd[5,5] = sum(beam.dist.z .* beam.dist.z)/beam.num_macro 
    beam.moment2nd[5,6] = sum(beam.dist.z .* beam.dist.dp)/beam.num_macro 

    beam.moment2nd[6,6] = sum(beam.dist.dp .* beam.dist.dp)/beam.num_macro 

    
    beam.moment2nd .= beam.moment2nd+beam.moment2nd'-diagm(diag(beam.moment2nd))

end

function get_emittance!(beam::BunchedBeam)
    # calculate emittance
    get_centroid!(beam)
    get_2nd_moments!(beam)
    @inbounds for i in 1:6
        for j in 1:6
            beam.moment2nd[i,j] = beam.moment2nd[i,j] - beam.centroid[i]*beam.centroid[j]
        end
    end
    
    @inbounds for i in 1:3
        beam.emittance[i] = sqrt(det(@view beam.moment2nd[2*i-1:2*i,2*i-1:2*i]))
    end
    beam.emittance[3] = beam.emittance[3] * beam.beta * beam.total_energy / 2.99792458e8

end

function initilize_6DGaussiandist!(beam::BunchedBeam, optics::AbstractOptics4D, rf::AbstractCavity, αc::Float64)
    # 6D Gaussian distribution
    beam.dist.x .= randn(beam.num_macro)
    beam.dist.px .= randn(beam.num_macro)
    beam.dist.y .= randn(beam.num_macro)
    beam.dist.py .= randn(beam.num_macro)
    beam.dist.z .= randn(beam.num_macro)
    beam.dist.dp .= randn(beam.num_macro)

    get_centroid!(beam)

    beam.dist.x .= beam.dist.x .- beam.centroid[1]
    beam.dist.px .= beam.dist.px .- beam.centroid[2]
    beam.dist.y .= beam.dist.y .- beam.centroid[3]
    beam.dist.py .= beam.dist.py .- beam.centroid[4]
    beam.dist.z .= beam.dist.z .- beam.centroid[5]
    beam.dist.dp .= beam.dist.dp .- beam.centroid[6]

    get_2nd_moments!(beam)
    
    eigval,eigvec=eigen(beam.moment2nd)
    
    mscale=eigvec * diagm(1.0 ./ sqrt.(eigval)) * eigvec'
    # # #beam.dist=mscale*beam.dist
    newx = mscale[1,1] .* beam.dist.x + mscale[1,2] .* beam.dist.px + mscale[1,3] .* beam.dist.y + mscale[1,4] .* beam.dist.py + mscale[1,5] .* beam.dist.z + mscale[1,6] .* beam.dist.dp
    newpx = mscale[2,1] .* beam.dist.x + mscale[2,2] .* beam.dist.px + mscale[2,3] .* beam.dist.y + mscale[2,4] .* beam.dist.py + mscale[2,5] .* beam.dist.z + mscale[2,6] .* beam.dist.dp
    newy = mscale[3,1] .* beam.dist.x + mscale[3,2] .* beam.dist.px + mscale[3,3] .* beam.dist.y + mscale[3,4] .* beam.dist.py + mscale[3,5] .* beam.dist.z + mscale[3,6] .* beam.dist.dp
    newpy = mscale[4,1] .* beam.dist.x + mscale[4,2] .* beam.dist.px + mscale[4,3] .* beam.dist.y + mscale[4,4] .* beam.dist.py + mscale[4,5] .* beam.dist.z + mscale[4,6] .* beam.dist.dp
    newz = mscale[5,1] .* beam.dist.x + mscale[5,2] .* beam.dist.px + mscale[5,3] .* beam.dist.y + mscale[5,4] .* beam.dist.py + mscale[5,5] .* beam.dist.z + mscale[5,6] .* beam.dist.dp
    newdp = mscale[6,1] .* beam.dist.x + mscale[6,2] .* beam.dist.px + mscale[6,3] .* beam.dist.y + mscale[6,4] .* beam.dist.py + mscale[6,5] .* beam.dist.z + mscale[6,6] .* beam.dist.dp
    
    
    beam.dist.x .= newx .* sqrt(beam.emittance[1]*optics.optics_x.beta)
    beam.dist.px .= newpx .* sqrt(beam.emittance[1]/optics.optics_x.beta)
    beam.dist.y .= newy .* sqrt(beam.emittance[2]*optics.optics_y.beta)
    beam.dist.py .= newpy .* sqrt(beam.emittance[2]/optics.optics_y.beta)
    beam.dist.px .+= beam.dist.x .* (optics.optics_x.alpha/optics.optics_x.beta)
    beam.dist.py .+= beam.dist.y .* (optics.optics_y.alpha/optics.optics_y.beta)

    
    # #generate longtiudinal distribution based on small amplitude approximation
    eta_p=αc-1.0/beam.gamma^2
    
    Qs=sqrt(rf.v*rf.h*abs(eta_p*cos(rf.ϕs))/2/π/beam.beta^2/beam.total_energy)

    emit_deltap_z=beam.emittance[3]*2.99792458e8/beam.beta/beam.total_energy
    invbeta_deltap_z=Qs*rf.k/rf.h/abs(eta_p)
    beam.dist.z .= newz .*sqrt(emit_deltap_z/invbeta_deltap_z)
    beam.dist.dp .= newdp .*sqrt(emit_deltap_z*invbeta_deltap_z)
    
    
end




function histogram1DinZ!(beam::BunchedBeam, nbins::Int64)
    # histogram in z
    zhist=zeros(nbins)
    zhist_edges=zeros(nbins+1)
    zmax=maximum(beam.dist.z)
    zmin=minimum(beam.dist.z)
    zhist_edges .= collect(range(zmin-(zmax-zmin)/nbins, zmax+(zmax-zmin)/nbins, length=nbins+1))
    zsep=(zhist_edges[end]-zhist_edges[1])/nbins
    @inbounds for i in 1:beam.num_macro
        ibin = (beam.dist.z[i] - zhist_edges[1])/zsep  # number of bin from 0
        dx = round(ibin) - ibin
        binnum=Int64(floor(ibin)+1)
        beam.inzindex[i] = binnum 
        neighbor = binnum + Int64(sign(dx))
        ratio = (0.5 - abs(dx))/0.5
        weight_neighbor = 0.5 * ratio^2
        zhist[binnum] += 1.0 - weight_neighbor
        zhist[neighbor] += weight_neighbor
    end
    return zhist, zhist_edges
end


        




