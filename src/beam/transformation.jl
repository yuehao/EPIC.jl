using StaticArrays


abstract type AbstractTransformation end
abstract type AbstractTransformation2D <: AbstractTransformation end
abstract type AbstractTransformation4D <: AbstractTransformation end


struct LinearTransformation2D <: AbstractTransformation2D
    beta::Float64
    alpha::Float64
    gamma::Float64
    phaseadv::Float64
    eta::Float64
    etap::Float64
    umap::SMatrix{2,2}
    invumap::SMatrix{2,2}
    Optics2D(b::Float64, a::Float64, d::Float64, dp::Float64)=new(b,a,(1+a*a)/b,0,d,dp, 
                                                                @SMatrix [1/sqrt(beta) 0 ; alpha/sqrt(beta)  sqrt(beta)],
                                                                @SMatrix [sqrt(beta) 0 ; -alpha*sqrt(beta)  1/sqrt(beta)])
end

struct LinearTransformation4D <: AbstractTransformation4D   # 4D linear transformation Uncoupled
    linearmap_x::LinearTransformation2D
    linearmap_y::LinearTransformation2D
    umap::SMatrix{4,4}
    invumap::SMatrix{4,4}
    Optics4D(lmx::LinearTransformation2D, lmy::LinearTransformation2D)=new(lmx,lmy,
                                                                @SMatrix [lmx.umap zeros(2,2); zeros(2,2) lmy.umap],
                                                                @SMatrix [lmx.invumap zeros(2,2); zeros(2,2) lmy.invumap])
end

struct TransCouplingTransformation <: AbstractTransformation4D
    uncoupled4D::LinearTransformation4D
    coupling::SMatrix{4,4}
    Optics4D(uncoupled4D::LinearTransformation4D, coupling::SMatrix{4,4})=new(uncoupled4D,coupling)
end


struct Optics2DUncoupled <: AbstractOptics2D
    od1::Optics1D
    od2::Optics1D
    Optics2DUncoupled(o1::Optics1D, o2::Optics1D) = new(o1,o2)
end

struct Optics2DEdTeng <: AbstractOptics2D
    od1::Optics1D
    od2::Optics1D
    α::Float64
    Rotm::SMatrix{2, 2, Float64}
end



function normalization(o2::Optics2DUncoupled)
    m1=normalization(o2.od1)
    m2=normalization(o2.od2)
    zm=@SMatrix [0.0 0.0; 0.0 0.0]

    return SMatrix{4,4}([m1 zm; zm m2])
end


function invnormalization(o2::Optics2DUncoupled)
    m1=invnormalization(o2.od1)
    m2=invnormalization(o2.od2)
    zm=@SMatrix [0.0 0.0; 0.0 0.0]
    return SMatrix{4,4}([m1 zm; zm m2])
end



end