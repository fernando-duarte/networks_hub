import Pkg; Pkg.activate("optim"); Pkg.instantiate()
#Pkg.resolve();Pkg.gc();Pkg.precompile();
#include("NetworkType.jl"); using .NetworkType

using DataFrames
include("NetDefs.jl"); using .NetDefs

include("NonLinProbPrecompile.jl"); 
using .PrecompileExtpow, .NonLinProbPrecompileSumpContraction, .NonLinProbPrecompileContractionQuadUnif, .NonLinProbPrecompileNum, .NonLinProbPrecompileObj

include("ExtendedPowerDist.jl"); #include("IncBetaDer.jl"); #using .IncBetaDer

using ModelingToolkit, Symbolics, SparseArrays, LinearAlgebra
