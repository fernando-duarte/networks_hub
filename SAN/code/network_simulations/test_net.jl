import Pkg;Pkg.activate("nlp")
Pkg.instantiate()
using LinearAlgebra, DataFrames, XLSX, Missings, Test, SpecialFunctions,  SparseArrays , Random, BenchmarkTools, Distributions, JLD2
using JuMP,Ipopt, NLsolve, Optim, LineSearches
using NLPModels, ADNLPModels, Percival, Nonconvex, Hyperopt
using ForwardDiff
using BSON, DataFrames, NLsolve, Test
include("IncBetaDer.jl")
include("NetOptimProb.jl")

## load data

# xf = XLSX.readxlsx("glasserman_young_example.xlsx") 
# data = vcat( [(XLSX.eachtablerow(xf[s]) |> DataFrames.DataFrame) for s in XLSX.sheetnames(xf)]... ) #for s in XLSX.sheetnames(xf) if (s!="Aggregates Composition" && s!="Dealer Aggregates" && s!="Approx Aggregates")
# unique!(data) # delete duplicate rows, use `nonunique(data)` to see if there are any duplicates
# data = data[isequal.(data.qt_dt,195), :] # keep quarter == 195 = 2008q4
# sort!(data, :assets, rev = true)
# units = 1;
# data[:,[:w, :c, :assets, :p_bar, :b]] .= data[!,[:w, :c, :assets, :p_bar, :b]]./units
# # data.b[:] .= missing
# # data.c[:] .= missing

# names(data) # column names
# describe(data)
# show(data, allcols = true)
# bson("gy.bson",Dict(:data=>data))

# data_dict = BSON.load("gy.bson")
# show(data_dict[:data], allcols = true)


@testset "Example from Glasserman and Young" begin
    data_dict = BSON.load("gy.bson")
    N=5
    T = Float64
    b = [deepcopy(Array{T}(data_dict[:data][1:N,[:b]]))...]
    c = [deepcopy(Array{T}(data_dict[:data][1:N,[:c]]))...]
    w = [deepcopy(Array{T}(data_dict[:data][1:N,[:w]]))...]
    assets= [deepcopy(Array{T}(data_dict[:data][1:N,[:assets]]))...]
    p_bar = [deepcopy(Array{T}(data_dict[:data][1:N,[:p_bar]]))...]

    γ = T(0)
    mini_batch = 1;
    M = N^2+4*N+mini_batch*N
    
    y = 0; # parameter controlling degree to network linkages
    A  = [0 y/p_bar[1] 0 0 0; 0 0 0 y/p_bar[1] 0; 10.0/p_bar[3] 10.0/p_bar[3] 0 10.0/p_bar[3] 10.0/p_bar[3]; 0 0 0 0 y/p_bar[1]; y/p_bar[1] 0 0 0 0]; # matrix of relative liabilities 
    x = [0.0, 0.0, 94.0, 0.0, 0.0] # shock to the central node of size 94
    
    function clearing_vec!(F, p, x)
        F .= p.-min.(p_bar,max.(A'*p .+ c .- x,0))
    end
    sol = nlsolve((F, p)->clearing_vec!(F, p, x),p_bar, autodiff = :forward) 
    p= sol.zero; # the solution to the system
    @test converged(sol)

    sf = sum(p_bar-p) # payment shortfalls
    dl = sum(x) # direct losses
    L = sf+dl

    @test sf==88
    @test dl==94
    @test L==182

    sf_d = sum(max.(0,x-w)) # payment shortfalls for disconnected network
    dl_d = sum(x) # direct losses for disconnected network
    L_d = sf_d+dl_d # expected spillovers for disconnected network

    @test sf_d==84
    @test dl_d==94
    @test L_d==178

    α = ones(N) # arbitrary value, plays no role in this one-shock example
    β = ones(N) # arbitrary value, plays no role in this one-shock example
    z = vcat(A...,b...,c...,α,β,p_bar)
    out = NetOptimProb.obj_full(z,x./c)
    @test L== NetOptimProb. obj_full(z,x./c)
    @test L_d== NetOptimProb.loss_d(z,x./c)   
end

