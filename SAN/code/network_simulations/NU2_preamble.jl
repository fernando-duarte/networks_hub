# module NetworkUtils2

# using Convex, SCS, NLsolve
# using LinearAlgebra, Distributions, Random, SpecialFunctions,  SparseArrays, MeasureTheory
# using BenchmarkTools, Test
# using DataFrames, BSON, XLSX, JLD2, Missings
# using Dates, MonthlyDates
# # using Plots, Gadfly, Cairo, Fontconfig, LightGraphs, SimpleWeightedGraphs, GraphRecipes
# using Setfield, Parameters, Suppressor, Flatten
# using MathOptInterface
# const MOI = MathOptInterface
# using Reexport

# include("IncBetaDer.jl")
# using .IncBetaDer
# @reexport using .NetworkType
# export optimParam, c_lin, c_quad, c_chance, c_p, c_pNL, c_, c_NL, c_, c_NL, obj, obj_full, loss_d

import Pkg; Pkg.activate("optim"); Pkg.instantiate()
#Pkg.resolve();Pkg.gc();Pkg.precompile();
using BSON, DataFrames, LinearAlgebra, SparseArrays, Parameters, Statistics, NLsolve, Quadrature, Cuba,  DistributionsAD, SpecialFunctions, BlockArrays, Revise, Test
using ModelingToolkit, GalacticOptim, Optim, NonlinearSolve, Symbolics
import Distributions
include("NetworkType.jl"); using .NetworkType
#include("IncBetaDer.jl"); #using .IncBetaDer
include("NetDefs.jl"); using .NetDefs
include("NonLinProbPrecompile.jl"); 
using .PrecompileExtpow, .NonLinProbPrecompileSumpContraction, .NonLinProbPrecompileContractionQuadUnif, .NonLinProbPrecompileNum, .NonLinProbPrecompile
include("ExtendedPowerDist.jl");



