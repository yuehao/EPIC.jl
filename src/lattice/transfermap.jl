
struct TransferMap2D <: AbstractTransferMap
    dim::Int64
    linearmap::SMatrix{2,2,Float64}
    TransferMap2D(mat::AbstractMatrix)=new(2, mat)
end

function TransferMap2D(o1::AbstractOptics2D, o2::AbstractOptics2D, phi::Float64)
    s,c=sincos(phi)
    rotation=@SMatrix [c s; -s c]
    return TransferMap1D(invnormal_mat(o2)*rotation*normal_mat(o1))
end

function TransferMap2D(o1::AbstractOptics2D, phi::Float64)
    return TransferMap1D(o1, o1, phi)
end

struct TransferMap4D <: AbstractTransferMap
    dim::Int64
    linearmap::SMatrix{4,4,Float64}
    TransferMap4D(mat::AbstractMatrix)=new(4, mat)
end
function TransferMap4D(o1::AbstractOptics4D, o2::AbstractOptics4D, phi1::Float64, phi2::Float64)
    s1,c1=sincos(phi1)
    r1=@SMatrix [c1 s1; -s1 c1]
    s2,c2=sincos(phi2)
    r2=@SMatrix [c2 s2; -s2 c2]
    zm=@SMatrix [0.0 0.0; 0.0 0.0]
    rotation=SMatrix{4,4}([r1 zm; zm r2])
    return TransferMap4D(invnormal_mat(o2)*rotation*normal_mat(o1))
end

function TransferMap4D(o1::AbstractOptics4D, phi1::Float64, phi2::Float64)
    return TransferMap4D(o1, o1, phi1, phi2)
end

function track!(ps6dcoor::AbstractVector{ps6d{T}}, tm::TransferMap4D) where T
    newx=tm.linearmap[1,1] .* ps6dcoor.x + tm.linearmap[1,2] .* ps6dcoor.px + tm.linearmap[1,3] .* ps6dcoor.y + tm.linearmap[1,4] .* ps6dcoor.py 
    newpx=tm.linearmap[2,1] .* ps6dcoor.x + tm.linearmap[2,2] .* ps6dcoor.px + tm.linearmap[2,3] .* ps6dcoor.y + tm.linearmap[2,4] .* ps6dcoor.py
    newy=tm.linearmap[3,1] .* ps6dcoor.x + tm.linearmap[3,2] .* ps6dcoor.px + tm.linearmap[3,3] .* ps6dcoor.y + tm.linearmap[3,4] .* ps6dcoor.py
    ps6dcoor.py .= tm.linearmap[4,1] .* ps6dcoor.x + tm.linearmap[4,2] .* ps6dcoor.px + tm.linearmap[4,3] .* ps6dcoor.y + tm.linearmap[4,4] .* ps6dcoor.py
    ps6dcoor.x .= newx
    ps6dcoor.px .= newpx
    ps6dcoor.y .= newy
    
end


function track!(beam::BunchedBeam, tm::TransferMap4D)
    @inbounds for i in 1:beam.num_macro
        beam.dist.x[i], beam.dist.px[i], beam.dist.y[i], beam.dist.py[i] = tm.linearmap * [beam.dist.x[i], beam.dist.px[i], beam.dist.y[i], beam.dist.py[i]]
    end
end



