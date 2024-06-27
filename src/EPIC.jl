module EPIC
using StaticArrays
using StructArrays

include("abstypes.jl")
include("beam/beam.jl")

include("lattice/elements.jl")
include("lattice/lattice.jl")
include("optics/optics.jl")

export PROTON, ELECTRON, GOLDION
export ps6d

export BunchedBeam
export get_centroid!,get_2nd_moment!, get_emittance!
export initilize_6DGaussiandist!, histogram1DinZ!

export optics4DUC

export StrongGaussianBeam, StrongThinGaussianBeam, initilize_zslice!, Bassetti_Erskine!, crab_crossing_setup!
export SpaceChargeLens
export Drift, ThinCorrector
export TransferMap4D, TransferMap4DChrom, LongitudinalRFMap, OneTurnRadiation
export AccelCavity, CrabCavity, easyCrabCavity, LongitudinalRLCWake, LongitudinalWake
export LorentzBoost, InvLorentzBoost
export Lattice
export track!
export get_synchrotron_tune

end