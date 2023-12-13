abstract type AbstractBeam end


abstract type AbstractOptics end
abstract type AbstractOptics2D <:AbstractOptics end
abstract type AbstractOptics4D <:AbstractOptics end


abstract type AbstractElement end

abstract type AbstractCavity <:AbstractElement end

abstract type AbstractWakefield <:AbstractElement end
abstract type AbstractLongiWakefield <:AbstractWakefield end
abstract type AbstractTransWakefield <:AbstractWakefield end

abstract type AbstractRadiation <: AbstractElement end

abstract type AbstractTransferMap <:AbstractElement end

abstract type AbstractLattice end