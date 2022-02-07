struct TransferMap1D <: AbstractTransferMap
    dim::Int64
    linearmap::SMatrix{2,2,Float64}
    TransferMap1D(mat::AbstractMatrix)=new(2, mat)
end

function TransferMap1D(o1::Optics1D, o2::Optics1D, phi::Float64)
    s,c=sincos(phi)
    rotation=@SMatrix [c s; -s c]
    return TransferMap1D(invnormalization(o2)*rotation*normalization(o1))
end

function TransferMap1D(o1::Optics1D, phi::Float64)
    return TransferMap1D(o1, o1, phi)
end

struct TransferMap2D <: AbstractTransferMap
    dim::Int64
    linearmap::SMatrix{4,4,Float64}
    TransferMap2D(mat::AbstractMatrix)=new(4, mat)
end
function TransferMap2D(o1::Optics2DUncoupled, o2::Optics2DUncoupled, phi1::Float64, phi2::Float64)
    s1,c1=sincos(phi1)
    r1=@SMatrix [c1 s1; -s1 c1]
    s2,c2=sincos(phi2)
    r2=@SMatrix [c2 s2; -s2 c2]
    zm=@SMatrix [0.0 0.0; 0.0 0.0]
    rotation=SMatrix{4,4}([r1 zm; zm r2])
    return TransferMap2D(invnormalization(o2)*rotation*normalization(o1))
end

function TransferMap2D(o1::Optics2D, phi1::Float64, phi2::Float64)
    return TransferMap2D(o1, o1, phi1, phi2)
end

function track(coor6d::AbstractVector, tm::AbstractTransferMap)
    coorview=@view coor6d[1:tm.dim]
    coorview=tm.linearmap*coorview
    return nothing
end

