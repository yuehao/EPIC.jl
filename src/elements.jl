import .optics
using StaticArrays


module elements
abstract type AbstractElement end
abstract type AbstractCavity <:AbstractElement end

abstract type AbstractTransferMap <:AbstractElement end
abstract type AbstractStrongBeamBeam <:AbstractElement end

include("cavity.jl")
include("transfermap.jl")
include("strongbb.jl")


end 