using LinearAlgebra


const distribution1D_type=Set([:custom, :gaussian, :uniform])
const zslice_type=Set([:evenzsep, :evennpar])


function initilize_6DGaussiandist!(beam::BunchedBeam, optics::AbstractOptics4D, rf::AbstractCavity, αc::Float64)
    # 6D Gaussian distribution

    moment_matrix=zeros(6,6)

    moment_matrix[1,1] = sum(beam.dist.x .* beam.dist.x)/beam.num_macro
    moment_matrix[1,2] = sum(beam.dist.x .* beam.dist.px)/beam.num_macro
    moment_matrix[1,3] = sum(beam.dist.x .* beam.dist.y)/beam.num_macro
    moment_matrix[1,4] = sum(beam.dist.x .* beam.dist.py)/beam.num_macro
    moment_matrix[1,5] = sum(beam.dist.x .* beam.dist.z)/beam.num_macro
    moment_matrix[1,6] = sum(beam.dist.x .* beam.dist.dp)/beam.num_macro

    moment_matrix[2,2] = sum(beam.dist.px .* beam.dist.px)/beam.num_macro
    moment_matrix[2,3] = sum(beam.dist.px .* beam.dist.y)/beam.num_macro
    moment_matrix[2,4] = sum(beam.dist.px .* beam.dist.py)/beam.num_macro
    moment_matrix[2,5] = sum(beam.dist.px .* beam.dist.z)/beam.num_macro
    moment_matrix[2,6] = sum(beam.dist.px .* beam.dist.dp)/beam.num_macro

    moment_matrix[3,3] = sum(beam.dist.y .* beam.dist.y)/beam.num_macro
    moment_matrix[3,4] = sum(beam.dist.y .* beam.dist.py)/beam.num_macro
    moment_matrix[3,5] = sum(beam.dist.y .* beam.dist.z)/beam.num_macro
    moment_matrix[3,6] = sum(beam.dist.y .* beam.dist.dp)/beam.num_macro

    moment_matrix[4,4] = sum(beam.dist.py .* beam.dist.py)/beam.num_macro
    moment_matrix[4,5] = sum(beam.dist.py .* beam.dist.z)/beam.num_macro
    moment_matrix[4,6] = sum(beam.dist.py .* beam.dist.dp)/beam.num_macro

    moment_matrix[5,5] = sum(beam.dist.z .* beam.dist.z)/beam.num_macro
    moment_matrix[5,6] = sum(beam.dist.z .* beam.dist.dp)/beam.num_macro

    moment_matrix[6,6] = sum(beam.dist.dp .* beam.dist.dp)/beam.num_macro


    
    moment_matrix=moment_matrix+moment_matrix'-diagm(diag(moment_matrix))
    eigval,eigvec=eigen(moment_matrix)
    
    mscale=eigvec * diagm(1.0 ./ sqrt.(eigval)) * eigvec'
    # # #beam.dist=mscale*beam.dist
    newx = mscale[1,1] .* beam.dist.x + mscale[1,2] .* beam.dist.px + mscale[1,3] .* beam.dist.y + mscale[1,4] .* beam.dist.py + mscale[1,5] .* beam.dist.z + mscale[1,6] .* beam.dist.dp
    newpx = mscale[2,1] .* beam.dist.x + mscale[2,2] .* beam.dist.px + mscale[2,3] .* beam.dist.y + mscale[2,4] .* beam.dist.py + mscale[2,5] .* beam.dist.z + mscale[2,6] .* beam.dist.dp
    newy = mscale[3,1] .* beam.dist.x + mscale[3,2] .* beam.dist.px + mscale[3,3] .* beam.dist.y + mscale[3,4] .* beam.dist.py + mscale[3,5] .* beam.dist.z + mscale[3,6] .* beam.dist.dp
    newpy = mscale[4,1] .* beam.dist.x + mscale[4,2] .* beam.dist.px + mscale[4,3] .* beam.dist.y + mscale[4,4] .* beam.dist.py + mscale[4,5] .* beam.dist.z + mscale[4,6] .* beam.dist.dp
    newz = mscale[5,1] .* beam.dist.x + mscale[5,2] .* beam.dist.px + mscale[5,3] .* beam.dist.y + mscale[5,4] .* beam.dist.py + mscale[5,5] .* beam.dist.z + mscale[5,6] .* beam.dist.dp
    newdp = mscale[6,1] .* beam.dist.x + mscale[6,2] .* beam.dist.px + mscale[6,3] .* beam.dist.y + mscale[6,4] .* beam.dist.py + mscale[6,5] .* beam.dist.z + mscale[6,6] .* beam.dist.dp
    
    beam.dist.x .= newx
    beam.dist.px .= newpx
    beam.dist.y .= newy
    beam.dist.py .= newpy
    beam.dist.z .= newz
    beam.dist.dp .= newdp

    beam.dist.x .= beam.dist.x .* sqrt(beam.emittance[1]*optics.optics_x.beta)
    beam.dist.px .= beam.dist.px .* sqrt(beam.emittance[1]/optics.optics_x.beta)
    beam.dist.px .+= beam.dist.x .* (optics.optics_x.alpha/optics.optics_x.beta)

    beam.dist.y .= beam.dist.y .* sqrt(beam.emittance[2]*optics.optics_y.beta)
    beam.dist.py .= beam.dist.py .* sqrt(beam.emittance[2]/optics.optics_y.beta)
    beam.dist.py .+= beam.dist.y .* (optics.optics_y.alpha/optics.optics_y.beta)

    # #generate longtiudinal distribution based on small amplitude approximation
    eta_p=αc-1.0/beam.gamma^2
    
    Qs=sqrt(rf.v*rf.h*abs(eta_p*cos(rf.ϕs))/2/π/beam.beta^2/beam.total_energy)

    emit_deltap_z=beam.emittance[3]*2.99792458e8/beam.beta/beam.total_energy
    invbeta_deltap_z=Qs*rf.k/rf.h/abs(eta_p)
    beam.dist.z .= beam.dist.z .*sqrt(emit_deltap_z/invbeta_deltap_z)
    beam.dist.dp .= beam.dist.dp .*sqrt(emit_deltap_z*invbeta_deltap_z)
    
    
end
