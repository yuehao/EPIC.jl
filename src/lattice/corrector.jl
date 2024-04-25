struct ThinCorrector <: AbstractCorrector
    dx::Float64
    dy::Float64
end

function track!(ps6dcoor::AbstractVector{ps6d{T}}, tc::ThinCorrector) where T
    @inbounds Threads.@threads for i in eachindex(ps6dcoor.x)
        ps6dcoor.px[i] += tc.dx
        ps6dcoor.py[i] += tc.dy
    end
    nothing
end

function track_bcast!(ps6dcoor::AbstractVector{ps6d{T}}, tc::ThinCorrector) where T
    ps6dcoor.px .+= tc.dx
    ps6dcoor.py .+= tc.dy
    nothing
end

function track!(beam::BunchedBeam, tc::ThinCorrector)
    track!(beam.dist, tc)
    nothing
end