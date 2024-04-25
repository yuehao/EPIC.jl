
struct TransferMap2D <: AbstractTransverseMap
    dim::Symbol
    tune::Float64
    linearmap::SMatrix{2,2,Float64}
    TransferMap2D(direction::Symbol, tune::Float64, mat::AbstractMatrix) = new(direction, tune, mat)
end

function TransferMap2D(direction::Symbol, o1::AbstractOptics2D, o2::AbstractOptics2D, tune::Float64)
    if direction != :x && direction != :y
        error("direction should be :x or :y")
    end
    s,c=sincos(2π*tune)
    rotation=@SMatrix [c s; -s c]
    return TransferMap2D(direction, tune, invnormal_mat(o2)*rotation*normal_mat(o1))
end

function TransferMap2D(direction, o1::AbstractOptics2D, tune::Float64)
    return TransferMap2D(direction, o1, o1, tune)
end

function track!(ps6dcoor::AbstractVector{ps6d{T}}, tm::TransferMap2D, temp) where T
    if tm.dim==:x
        @inbounds Threads.@threads for i in eachindex(ps6dcoor.x)
            temp[i] = tm.linearmap[1,1] * ps6dcoor.x[i] + tm.linearmap[1,2] * ps6dcoor.px[i]
            ps6dcoor.px[i] = tm.linearmap[2,1] * ps6dcoor.x[i] + tm.linearmap[2,2] * ps6dcoor.px[i]
            ps6dcoor.x[i] = temp[i]
        end
        #temp  .= tm.linearmap[1,1] .* ps6dcoor.x .+ tm.linearmap[1,2] .* ps6dcoor.px
        #ps6dcoor.px .= tm.linearmap[2,1] .* ps6dcoor.x .+ tm.linearmap[2,2] .* ps6dcoor.px
        #ps6dcoor.x .= temp
    else
        @inbounds Threads.@threads for i in eachindex(ps6dcoor.y)
            temp[i] = tm.linearmap[1,1] * ps6dcoor.y[i] + tm.linearmap[1,2] * ps6dcoor.py[i]
            ps6dcoor.py[i] = tm.linearmap[2,1] * ps6dcoor.y[i] + tm.linearmap[2,2] * ps6dcoor.py[i]
            ps6dcoor.y[i] = temp[i]
        end
        #temp  .= tm.linearmap[1,1] .* ps6dcoor.y .+ tm.linearmap[1,2] .* ps6dcoor.py
        #ps6dcoor.py .= tm.linearmap[2,1] .* ps6dcoor.y .+ tm.linearmap[2,2] .* ps6dcoor.py
        #ps6dcoor.y .= temp
    end
    return nothing
end



function track!(beam::BunchedBeam, tm::TransferMap2D)
    trackx!(beam.dist, tm, beam.temp1)   
end

struct TransferMap2DChrom <: AbstractTransverseMap
    dim::Int64
    tune::Float64
    chrom::Float64
    umat::SMatrix{2,2,Float64, 4}
    invumat::SMatrix{2,2,Float64, 4}
    TransferMap2DChrom(tune, chrom, umat::AbstractMatrix, invumat::AbstractMatrix)=new(2, tune, chrom, umat, invumat)
end

function TransferMap2DChrom(o1::AbstractOptics2D, o2::AbstractOptics2D, tune::Float64, chrom::Float64)
    return TransferMap1DChrom(tune, chrom, normal_mat(o1), invnormal_mat(o2))
end

function TransferMap2DChrom(o1::AbstractOptics2D, tune::Float64, chrom::Float64)
    return TransferMap1DChrom(tune, chrom, normal_mat(o1), invnormal_mat(o1))
end



struct TransferMap4D <: AbstractTransverseMap
    dim::Int64
    tune::SVector{2,Float64}
    linearmap::SMatrix{4,4,Float64,16}
    TransferMap4D(tune::AbstractVector, mat::AbstractMatrix)=new(4, tune, mat)
end

function TransferMap4D(o1::AbstractOptics4D, o2::AbstractOptics4D, tunex::Float64, tuney::Float64)
    s1,c1=sincos(2π*tunex)
    r1=@SMatrix [c1 s1; -s1 c1]
    s2,c2=sincos(2π*tuney)
    r2=@SMatrix [c2 s2; -s2 c2]
    zm=@SMatrix [0.0 0.0; 0.0 0.0]
    rotation=SMatrix{4,4}([r1 zm; zm r2])
    tune=SVector{2,Float64}(tunex, tuney)
    return TransferMap4D(tune, invnormal_mat(o2)*rotation*normal_mat(o1))
end

function TransferMap4D(o1::AbstractOptics4D, tunex::Float64, tuney::Float64)
    return TransferMap4D(o1, o1, tunex, tuney)
end

function track!(ps6dcoor::AbstractVector{ps6d{T}}, tm::TransferMap4D, tempx, temppx, tempy) where T
    @inbounds Threads.@threads for i in eachindex(ps6dcoor.x)
        tempx[i] = tm.linearmap[1,1] * ps6dcoor.x[i] + tm.linearmap[1,2] * ps6dcoor.px[i] + tm.linearmap[1,3] * ps6dcoor.y[i] + tm.linearmap[1,4] * ps6dcoor.py[i]
        temppx[i] = tm.linearmap[2,1] * ps6dcoor.x[i] + tm.linearmap[2,2] * ps6dcoor.px[i] + tm.linearmap[2,3] * ps6dcoor.y[i] + tm.linearmap[2,4] * ps6dcoor.py[i]
        tempy[i] = tm.linearmap[3,1] * ps6dcoor.x[i] + tm.linearmap[3,2] * ps6dcoor.px[i] + tm.linearmap[3,3] * ps6dcoor.y[i] + tm.linearmap[3,4] * ps6dcoor.py[i]
        ps6dcoor.py[i] = tm.linearmap[4,1] * ps6dcoor.x[i] + tm.linearmap[4,2] * ps6dcoor.px[i] + tm.linearmap[4,3] * ps6dcoor.y[i] + tm.linearmap[4,4] * ps6dcoor.py[i]
        ps6dcoor.x[i] = tempx[i]
        ps6dcoor.px[i] = temppx[i]
        ps6dcoor.y[i] = tempy[i]
    end

    # tempx  .= tm.linearmap[1,1] .* ps6dcoor.x .+ tm.linearmap[1,2] .* ps6dcoor.px .+ tm.linearmap[1,3] .* ps6dcoor.y .+ tm.linearmap[1,4] .* ps6dcoor.py 
    # temppx .= tm.linearmap[2,1] .* ps6dcoor.x .+ tm.linearmap[2,2] .* ps6dcoor.px .+ tm.linearmap[2,3] .* ps6dcoor.y .+ tm.linearmap[2,4] .* ps6dcoor.py
    # tempy  .= tm.linearmap[3,1] .* ps6dcoor.x .+ tm.linearmap[3,2] .* ps6dcoor.px .+ tm.linearmap[3,3] .* ps6dcoor.y .+ tm.linearmap[3,4] .* ps6dcoor.py
    # ps6dcoor.py .= tm.linearmap[4,1] .* ps6dcoor.x .+ tm.linearmap[4,2] .* ps6dcoor.px .+ tm.linearmap[4,3] .* ps6dcoor.y .+ tm.linearmap[4,4] .* ps6dcoor.py
    # ps6dcoor.x .= tempx
    # ps6dcoor.px .= temppx
    # ps6dcoor.y .= tempy
    
    return nothing
end


function track!(beam::BunchedBeam, tm::TransferMap4D)
    track!(beam.dist, tm,  beam.temp1, beam.temp2, beam.temp3)
end

struct TransferMap4DChrom <: AbstractTransverseMap
    dim::Int64
    tune::SVector{2,Float64}
    chrom::SVector{2,Float64}
    umat::SMatrix{4,4,Float64, 16}
    invumat::SMatrix{4,4,Float64, 16}
    TransferMap4DChrom(tune::AbstractVector, chrom::AbstractVector, umat::AbstractMatrix, invumat::AbstractMatrix)=new(4, tune, chrom, umat, invumat)
end
function TransferMap4DChrom(o1::AbstractOptics4D, o2::AbstractOptics4D, tunex::Float64, tuney::Float64, chromx::Float64, chromy::Float64)
    tune=SVector{2,Float64}(tunex, tuney)
    chrom=SVector{2,Float64}(chromx, chromy)
    return TransferMap4DChrom(tune, chrom, normal_mat(o1), invnormal_mat(o2))
end

function TransferMap4DChrom(o1::AbstractOptics4D, tunex::Float64, tuney::Float64, chromx::Float64, chromy::Float64)
    return TransferMap4DChrom(o1, o1, tunex, tuney, chromx, chromy)
end




function track!(ps6dcoor::AbstractVector{ps6d{T}}, tm::TransferMap4DChrom, temp1, temp2, temp3, sinphi, cosphi) where T

    @inbounds Threads.@threads for i in eachindex(temp1)
        temp1[i] = tm.umat[1,1] * ps6dcoor.x[i] + tm.umat[1,2] * ps6dcoor.px[i] + tm.umat[1,3] * ps6dcoor.y[i] + tm.umat[1,4] * ps6dcoor.py[i]
        temp2[i] = tm.umat[2,1] * ps6dcoor.x[i] + tm.umat[2,2] * ps6dcoor.px[i] + tm.umat[2,3] * ps6dcoor.y[i] + tm.umat[2,4] * ps6dcoor.py[i]
        temp3[i] = tm.umat[3,1] * ps6dcoor.x[i] + tm.umat[3,2] * ps6dcoor.px[i] + tm.umat[3,3] * ps6dcoor.y[i] + tm.umat[3,4] * ps6dcoor.py[i]
        ps6dcoor.py[i] = tm.umat[4,1] * ps6dcoor.x[i] + tm.umat[4,2] * ps6dcoor.px[i] + tm.umat[4,3] * ps6dcoor.y[i] + tm.umat[4,4] * ps6dcoor.py[i]
        ps6dcoor.x[i] = temp1[i]
        ps6dcoor.px[i] = temp2[i]
        ps6dcoor.y[i] = temp3[i]
    

        temp1[i] = 2π*tm.tune[1] + (2π*tm.chrom[1]) * ps6dcoor.dp[i]
        sinphi[i] = sin(temp1[i])
        cosphi[i] = cos(temp1[i])
        temp3[i] = cosphi[i] * ps6dcoor.x[i] + sinphi[i] * ps6dcoor.px[i]
        ps6dcoor.px[i] =  cosphi[i] * ps6dcoor.px[i] - sinphi[i] * ps6dcoor.x[i]
        ps6dcoor.x[i] = temp3[i]

        temp2[i] = 2π*tm.tune[2] + (2π*tm.chrom[2]) * ps6dcoor.dp[i]
        sinphi[i] = sin(temp2[i])
        cosphi[i] = cos(temp2[i])
        temp3[i] = cosphi[i] * ps6dcoor.y[i] + sinphi[i] * ps6dcoor.py[i]
        ps6dcoor.py[i] = cosphi[i] * ps6dcoor.py[i] - sinphi[i] * ps6dcoor.y[i]
        ps6dcoor.y[i] = temp3[i]

        temp1[i] = tm.invumat[1,1] * ps6dcoor.x[i] + tm.invumat[1,2] * ps6dcoor.px[i] + tm.invumat[1,3] * ps6dcoor.y[i] + tm.invumat[1,4] * ps6dcoor.py[i]
        temp2[i] = tm.invumat[2,1] * ps6dcoor.x[i] + tm.invumat[2,2] * ps6dcoor.px[i] + tm.invumat[2,3] * ps6dcoor.y[i] + tm.invumat[2,4] * ps6dcoor.py[i]
        temp3[i] = tm.invumat[3,1] * ps6dcoor.x[i] + tm.invumat[3,2] * ps6dcoor.px[i] + tm.invumat[3,3] * ps6dcoor.y[i] + tm.invumat[3,4] * ps6dcoor.py[i]
        ps6dcoor.py[i] = tm.invumat[4,1] * ps6dcoor.x[i] + tm.invumat[4,2] * ps6dcoor.px[i] + tm.invumat[4,3] * ps6dcoor.y[i] + tm.invumat[4,4] * ps6dcoor.py[i]
        ps6dcoor.x[i] = temp1[i]
        ps6dcoor.px[i] = temp2[i]
        ps6dcoor.y[i] = temp3[i]
    end

    # temp1 .= tm.umat[1,1] .* ps6dcoor.x .+ tm.umat[1,2] .* ps6dcoor.px .+ tm.umat[1,3] .* ps6dcoor.y .+ tm.umat[1,4] .* ps6dcoor.py 
    # temp2 .= tm.umat[2,1] .* ps6dcoor.x .+ tm.umat[2,2] .* ps6dcoor.px .+ tm.umat[2,3] .* ps6dcoor.y .+ tm.umat[2,4] .* ps6dcoor.py
    # temp3 .= tm.umat[3,1] .* ps6dcoor.x .+ tm.umat[3,2] .* ps6dcoor.px .+ tm.umat[3,3] .* ps6dcoor.y .+ tm.umat[3,4] .* ps6dcoor.py
    # ps6dcoor.py .= tm.umat[4,1] .* ps6dcoor.x .+ tm.umat[4,2] .* ps6dcoor.px .+ tm.umat[4,3] .* ps6dcoor.y .+ tm.umat[4,4] .* ps6dcoor.py
    # ps6dcoor.x .= temp1
    # ps6dcoor.px .= temp2
    # ps6dcoor.y .= temp3

    # temp1 .= 2π*tm.tune[1] .+ (2π*tm.chrom[1]) .* ps6dcoor.dp
    # sinphi .= sin.(temp1)
    # cosphi .= cos.(temp1)
    # temp3  .= cosphi .* ps6dcoor.x .+ sinphi .* ps6dcoor.px
    # ps6dcoor.px .=  cosphi .* ps6dcoor.px .- sinphi .* ps6dcoor.x
    # ps6dcoor.x .= temp3

    # temp2 .= 2π*tm.tune[2] .+ (2π*tm.chrom[2]) .* ps6dcoor.dp
    # sinphi .= sin.(temp2)
    # cosphi .= cos.(temp2)
    # temp3  .= cosphi .* ps6dcoor.y .+ sinphi .* ps6dcoor.py
    # ps6dcoor.py .= cosphi .* ps6dcoor.py .- sinphi .* ps6dcoor.y
    # ps6dcoor.y .= temp3

    # temp1  .= tm.invumat[1,1] .* ps6dcoor.x .+ tm.invumat[1,2] .* ps6dcoor.px .+ tm.invumat[1,3] .* ps6dcoor.y .+ tm.invumat[1,4] .* ps6dcoor.py
    # temp2 .= tm.invumat[2,1] .* ps6dcoor.x .+ tm.invumat[2,2] .* ps6dcoor.px .+ tm.invumat[2,3] .* ps6dcoor.y .+ tm.invumat[2,4] .* ps6dcoor.py
    # temp3  .= tm.invumat[3,1] .* ps6dcoor.x .+ tm.invumat[3,2] .* ps6dcoor.px .+ tm.invumat[3,3] .* ps6dcoor.y .+ tm.invumat[3,4] .* ps6dcoor.py
    # ps6dcoor.py .= tm.invumat[4,1] .* ps6dcoor.x .+ tm.invumat[4,2] .* ps6dcoor.px .+ tm.invumat[4,3] .* ps6dcoor.y .+ tm.invumat[4,4] .* ps6dcoor.py
    # ps6dcoor.x .= temp1
    # ps6dcoor.px .= temp2
    # ps6dcoor.y .= temp3

    return nothing

end

function track!(beam::BunchedBeam, tm::TransferMap4DChrom)
    track!(beam.dist, tm,  beam.temp1, beam.temp2, beam.temp3, beam.temp4, beam.temp5)
end



struct LongitudinalRFMap <: AbstractLongitudinalMap
    αc::Float64
    RF::AbstractAccelCavity
    LongitudinalRFMap(αc::Float64, RF::AbstractAccelCavity)=new(αc, RF)
end

function trackslippage!(ps6dcoor::AbstractVector{ps6d{T}}, alphac::Float64, h::Int64, k::Float64, gamma::Float64) where T
    eta = alphac - 1.0 / gamma / gamma
    # @inbounds Threads.@threads for i in eachindex(ps6dcoor.z)
    #     ps6dcoor.z[i] -= (2π * lm.RF.h * eta / lm.RF.k) * ps6dcoor.dp[i]
    # end
    ps6dcoor.z .-= (2π * h * eta / k) .* ps6dcoor.dp
    return nothing
end

function track!(beam::BunchedBeam, lm::LongitudinalRFMap)
    track!(beam, lm.RF)
    trackslippage!(beam.dist, lm.αc, lm.RF.h, lm.RF.k, beam.gamma)
end

function get_synchrotron_tune(beam::BunchedBeam, lmap::LongitudinalRFMap)
    eta = lmap.αc - 1.0 / beam.gamma / beam.gamma
    return sqrt(lmap.RF.v*lmap.RF.h*abs(eta*cos(lmap.RF.ϕs))/2/π/beam.beta^2/beam.total_energy)
    return sqrt(lm.RF.h * abs(lm.αc) / 2 / π / beam.beta / beam.total_energy)
end

struct LongitudinalCollectiveRFMap <: AbstractLongitudinalMap
    αc::Float64
    RF::AbstractAccelCavity
    wakefield::AbstractLongiWakefield
    nkick::Int64
    LongitudinalRFMap(αc::Float64, RF::AbstractAccelCavity, wake::AbstractLongiWakefield, nk::Int64)=new(αc, RF, wake, nk)
end

function track!(beam::BunchedBeam, lm::LongitudinalCollectiveRFMap)
    track!(beam, lm.RF)
    for i in 1:lm.nkick
        trackslippage!(beam.dist, lm.αc/lm.nkick, lm.RF.h, lm.RF.k, beam.gamma)
        histogram1DinZ!(beam)
        eN_b2E=beam.num_particle*1.6021766208e-19*beam.particle.charge^2/beam.total_energy/beam.beta/beam.beta/beam.particle.atomnum/m.nkick
        track!(beam.dist, rlcwake, beam.inzindex, eN_b2E, beam.znbin, beam.zhist, beam.zhist_edges, beam.ztemp1, beam.ztemp2, beam.ztemp3, beam.ztemp4)
    end
end