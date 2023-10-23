module BBMM

abstract type AbstractBeam end
include("beam/beam.jl")

abstract type AbstractOptics end


abstract type AbstractElement end
abstract type AbstractLattice end
include("lattice/lattice.jl")
include("lattice/elements.jl")





end