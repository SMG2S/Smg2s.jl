module Smg2s

using LinearAlgebra
using SparseArrays
using Random

export Nilpotent, Nilp
export Spectrum!, checkSpectrum
export initMat!
export nonherm, nonsym
export NilpxM, MxNilp

include("nilpotent.jl")
include("spectrum.jl")
include("init.jl")
include("nonsym.jl")
include("nonherm.jl")

end #endmodule
