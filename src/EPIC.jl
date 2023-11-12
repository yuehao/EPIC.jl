module EPIC
using StaticArrays
using StructArrays

include("abstypes.jl")
include("beam/beam.jl")

include("lattice/elements.jl")
include("optics/optics.jl")

export PROTON, ELECTRON, GOLDION

export BunchedBeam, initilize_6DGaussiandist!, histogram1DinZ!

export StrongGaussianBeam, initilize_zslice!
export optics4DUC, AccelCavity, LongitudinalRLCWake

export track!

end