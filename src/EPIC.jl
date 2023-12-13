module EPIC
using StaticArrays
using StructArrays

include("abstypes.jl")
include("beam/beam.jl")

include("lattice/elements.jl")
include("optics/optics.jl")

export PROTON, ELECTRON, GOLDION
export ps6d

export BunchedBeam
export get_centroid!,get_2nd_moments!, get_emittance!
export initilize_6DGaussiandist!, histogram1DinZ!

export optics4DUC

export StrongGaussianBeam, StrongThinGaussianBeam, initilize_zslice!, Bassetti_Erskine
export TransferMap4D, OneTurnRadiation
export AccelCavity, LongitudinalRLCWake

export track!

end