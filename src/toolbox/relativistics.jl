function getbeta(dpop::Float16,  βγ::Float64)
    p=(1+dpop)*βγ
    return p/sqrt(p*p+1)
end