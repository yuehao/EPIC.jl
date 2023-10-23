abstract type AbstractWakefield <:AbstractElement end

struct LongitudinalRLCWake <: AbstractWakefield
    freq::Float64
    Rshunt::Float64
    Q0::Float64
    LongitudinalRLCWake(freq,Rshunt,Q0)=new(freq,Rshunt,Q0)
end

function track!(beam::BunchedBeam, rlcwake::LongitudinalRLCWake)
    
end
