
using StaticArrays

abstract type AbstractElement end
abstract type AbstractCavity <:AbstractElement end
abstract type AbstractTransferMap <:AbstractElement end
abstract type AbstractBeamBeam <:AbstractElement end

abstract type AbstractWakefield <:AbstractElement end

include("cavity.jl")
include("transfermap.jl")
include("wakefield.jl")


