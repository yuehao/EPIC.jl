struct StrongBBThinGaussianBeam <: AbstractStrongBeamBeam
    amplitude::Float64
    rmssizex::Float64
    rmssizey::Float64
    zloc::Float64
    xoffset::Float64
    yoffset::Float64
    StrongBBThinGaussianBeam(amp::Float64, rx::Float64, ry::Float64, zloc::Float64=0.0, xoff::Float64=0.0, yoff::Float64=0.0)=new(amp,rx,ry,zloc,xoff,yoff)
end

struct StrongBBGaussianBeam <: AbstractStrongBeamBeam
    num_particle::Float64
    rmssizex::Float64
    rmssizey::Float64
    rmssizez::Float64
    opticsip::AbstractOptics2D
    xoffset::Float64
    yoffset::Float64

    StrongBBThinGaussianBeam(n::Float64, rx::Float64, ry::Float64, rz::Float64, optics::Opti
    xoff::Float64=0.0, yoff::Float64=0.0)=new(n,rx,ry,rz,xoff,yoff)
end


function Bassetti_Erskine(x::Float64, y::Float64, σx::Float64, σy::Float64)
	sqrtδσ2=sqrt(Complex(2*(σx*σx-σy*σy)))
	term1=erfcx(-1im*(x+1im*y)/sqrtδσ2)
	term2=erfcx(-1im*(x*σy/σx+1im*y*σx/σy)/sqrtδσ2)
	termexp=exp(-x*x/2/σx/σx-y*y/2/σy/σy)
	complex_e=-1im*2*sqrt(pi)/sqrtδσ2*(term1-termexp*term2)
	return real(complex_e), -imag(complex_e)
end

function track(coor6d::AbstractVector, stgb::StrongBBThinGaussianBeam)
    coor6d[1]+= (coor6d[2]*stgb.zloc)
    coor6d[3]+= (coor6d[4]*stgb.zloc)
    ex,ey=Bassetti_Erskine(coor6d[1]-stgb.xoffset, coor6d[3]-stgb.yoffset, stgb.rmssizex, stgb.rmssizey)
    coor6d[2]+=stgb.amplitude*ex
    coor6d[4]+=stgb.amplitude*ey
    coor6d[1]-= (coor6d[2]*stgb.zloc)
    coor6d[3]-= (coor6d[4]*stgb.zloc)
end

function track(coor6d::AbstractVector, stgb::StrongBBGaussianBeam)
    coor6d[1]+= (coor6d[2]*stgb.zloc)
    coor6d[3]+= (coor6d[4]*stgb.zloc)
    ex,ey=Bassetti_Erskine(coor6d[1]-stgb.xoffset, coor6d[3]-stgb.yoffset, stgb.rmssizex, stgb.rmssizey)
    coor6d[2]+=stgb.amplitude*ex
    coor6d[4]+=stgb.amplitude*ey
    coor6d[1]-= (coor6d[2]*stgb.zloc)
    coor6d[3]-= (coor6d[4]*stgb.zloc)
end


