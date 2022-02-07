using StaticArrays
module optics
export Optics1D, Optics2DUncoupled, Optics3DUncoupled, normalization, invnormalization
export AbstractOptics, AbstractOptics1D, AbstractOptics2D, AbstractOptics3D
abstract type AbstractOptics end
abstract type AbstractOptics1D <: AbstractOptics end
abstract type AbstractOptics2D <: AbstractOptics end
abstract type AbstractOptics3D <: AbstractOptics end

struct Optics1D <: AbstractOptics1D
    beta::Float64
    alpha::Float64
    gamma::Float64
    phaseadv::Float64
    eta::Float64
    etap::Float64
    Optics1D(b::Float64, a::Float64, d::Float64, dp::Float64)=new(b,a,(1+a*a)/b,0,d,dp)
end

struct Optics2DUncoupled <: AbstractOptics2D
    od1::Optics1D
    od2::Optics1D
    Optics2DUncoupled(o1::Optics1D, o2::Optics1D) = new(o1,o2)
end

struct Optics3DUncoupled <: AbstractOptics3D
    od1::Optics1D
    od2::Optics1D
    od3::Optics1D
    Optics2DUncoupled(o1::Optics1D, o2::Optics1D, o3::Optics1D) = new(o1,o2,o3)
end

function normalization(o1::Optics1D)
    sqrtbeta=sqrt(o1.beta)
    return @SMatrix [ 1.0/sqrtbeta 0 ; o1.alpha/sqrtbeta sqrtbeta]
end

function normalization(o2::Optics2DUncoupled)
    m1=normalization(o2.od1)
    m2=normalization(o2.od2)
    zm=@SMatrix [0.0 0.0; 0.0 0.0]

    return SMatrix{4,4}([m1 zm; zm m2])
end

function invnormalization(o1: Optics1D)
    sqrtbeta=sqrt(o1.beta)
    return @SMatrix [sqrtbeta 0 ; -o1.alpha*sqrtbetainv  sqrtbetainv]
end

function invnormalization(o2::Optics2DUncoupled)
    m1=invnormalization(o2.od1)
    m2=invnormalization(o2.od2)
    zm=@SMatrix [0.0 0.0; 0.0 0.0]

    return SMatrix{4,4}([m1 zm; zm m2])
end



end