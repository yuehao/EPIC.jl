
abstract type AbstractLongiWakefield <: AbstractWakefield end
abstract type AbstractTransWakefield <: AbstractWakefield end

struct LongitudinalRLCWake <: AbstractLongiWakefield
    freq::Float64
    Rshunt::Float64
    Q0::Float64
    wakefield::Function
end

function LongitudinalRLCWake(freq::Float64, Rshunt::Float64, Q0::Float64)
    Q0p=sqrt(Q0^2 - 1.0/4.0)
    ω0 = 2*pi*freq
    ω0p= ω0/Q0*Q0p
    wakefield = function (z::Float64)
        return Rshunt * ω0 /Q0 * (cos(ω0p * z) +  sin(ω0p * z) / 2 / Q0p) * exp(ω0 * z / 2 / Q0)
    end
    return LongitudinalRLCWake(freq, Rshunt, Q0, wakefield)
end

a=LongitudinalRLCWake(27e9, 50e3, 1.0)
zlist=collect(0.0: 1e-4: 1e-1)
wlist=map(a.wakefield, -zlist./2.99792458e8)
using Plots
plot(-zlist, wlist)

function track!(beam::BunchedBeam, rlcwake::LongitudinalRLCWake)
    
end
