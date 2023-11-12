
struct LongitudinalRLCWake <: AbstractLongiWakefield
    freq::Float64
    Rshunt::Float64
    Q0::Float64
    wakefield::Function
end

struct LongitudinalWake <: AbstractLongiWakefield
    wakefield::Function
end

struct TransverseWake <: AbstractTransWakefield
    wakefield::Function
end


function LongitudinalRLCWake(freq::Float64, Rshunt::Float64, Q0::Float64)
    Q0p=sqrt(Q0^2 - 1.0/4.0)
    ω0 = 2*pi*freq
    ω0p= ω0/Q0*Q0p
    wakefield = function (t::Float64)
        t>0 && return 0.0
        return Rshunt * ω0 /Q0 * (cos(ω0p * t) +  sin(ω0p * t) / 2 / Q0p) * exp(ω0 * t / 2 / Q0)
    end
    return LongitudinalRLCWake(freq, Rshunt, Q0, wakefield)
end

function track!(beam::BunchedBeam, rlcwake::AbstractLongiWakefield)
    nbins=Int64(round(sqrt(beam.num_macro)))÷2
    zhist, zhist_edges=histogram1DinZ!(beam, nbins)
    zhist_center=(zhist_edges[1:end-1]+zhist_edges[2:end])/2.0
    wakefield=rlcwake.wakefield.((zhist_center.-zhist_center[end]) ./ 2.99792458e8)
    
    wakepotential=zeros(nbins)
    eN_b2E=beam.num_particle*1.6021766208e-19*beam.particle.charge^2/beam.total_energy/beam.beta/beam.beta/beam.particle.atomnum
    @inbounds for i=1:nbins
        for j=i:nbins
            wakepotential[i]+=zhist[j]*wakefield[nbins-j+i]/beam.num_macro
        end
    end
    # wakepotential=DSP.conv(zhist, wakefield)[nbins:end]
    wakeatedge=zhist_edges.*0.0
    wakeatedge[2:end-1]=(wakepotential[1:end-1]+wakepotential[2:end])/2.0
    wakeatedge[1]= 2*wakeatedge[2]-wakeatedge[3]
    wakeatedge[end]= 2*wakeatedge[end-1]-wakeatedge[end-2]

    zsep=(zhist_edges[2]-zhist_edges[1])
    @inbounds for i in 1:beam.num_macro
        zloc=beam.dist.z[i]
        zindex=beam.inzindex[i]
        wake1=wakeatedge[zindex]
        wake2=wakeatedge[zindex+1]
        wakezloc=wake1+(wake2-wake1)*(zloc-zhist_edges[zindex])/zsep
        beam.dist.dp[i]-=wakezloc*eN_b2E
    end
end


