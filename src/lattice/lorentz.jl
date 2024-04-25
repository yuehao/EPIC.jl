

struct LorentzBoost <: AbstractLorentzBoost
    angle::Float64
    cosang::Float64
    tanang::Float64
    LorentzBoost(angle)=new(angle, cos(angle), tan(angle))
end

function track!(ps6dcoor::AbstractVector{ps6d{T}}, lb::LorentzBoost, mode::Symbol=:simple) where T
    if mode==:simple
        invcosang=1.0/lb.cosang
        @inbounds Threads.@threads for i in eachindex(ps6dcoor.x)
            ps6dcoor.x[i] += lb.tanang * ps6dcoor.z[i]
            ps6dcoor.dp[i] -= lb.tanang * ps6dcoor.px[i]
            ps6dcoor.px[i] *= 1.0/lb.cosang
            ps6dcoor.py[i] *= 1.0/lb.cosang
            ps6dcoor.z[i] *= 1.0/lb.cosang
        end
        
        #ps6dcoor.x .+= lb.tanang .* ps6dcoor.z
        #ps6dcoor.dp .-= lb.tanang .* ps6dcoor.px
        #ps6dcoor.px .*= invcosang
        #ps6dcoor.py .*= invcosang
        #ps6dcoor.z .*= invcosang
        return nothing
    end
    return nothing
end
function track!(beam::AbstractBeam, lb::LorentzBoost)
    track!(beam.dist, lb)
    return nothing
end




struct InvLorentzBoost <: AbstractLorentzBoost
    angle::Float64
    sinang::Float64
    cosang::Float64
    InvLorentzBoost(angle)=new(angle, sin(angle), cos(angle))
end

function track!(ps6dcoor::AbstractVector{ps6d{T}}, ilb::InvLorentzBoost, mode::Symbol=:simple) where T
    if mode==:simple
        @inbounds Threads.@threads for i in eachindex(ps6dcoor.x)
            ps6dcoor.x[i] -= ilb.sinang * ps6dcoor.z[i]
            ps6dcoor.dp[i] += ilb.sinang * ps6dcoor.px[i]
            ps6dcoor.px[i] *= ilb.cosang
            ps6dcoor.py[i] *= ilb.cosang
            ps6dcoor.z[i] *= ilb.cosang
        end
        #ps6dcoor.x .-= ilb.sinang .* ps6dcoor.z
        #ps6dcoor.dp .+= ilb.sinang .* ps6dcoor.px
        #ps6dcoor.px .*= ilb.cosang
        #ps6dcoor.py .*= ilb.cosang
        #ps6dcoor.z .*= ilb.cosang
        return nothing
    end
    return nothing
end
function track!(beam::AbstractBeam, ilb::InvLorentzBoost)
    track!(beam.dist, ilb)
    return nothing
end