## load packages
using Pkg
Pkg.add("FFTW")
Pkg.build("FFTW")
ENV["GRDIR"]="" ; Pkg.build("GR")
Pkg.add("PATHSolver")
Pkg.build("PATHSolver")
Pkg.add("SLEEFPirates")
Pkg.build("SLEEFPirates")
Pkg.add("HDF5")
Pkg.build("HDF5")

Pkg.add("DataFrames")
Pkg.add("XLSX")
Pkg.add("Missings")
Pkg.add("NLsolve")
Pkg.add("Complementarity")
Pkg.add("NLopt")
Pkg.add("JuMP")
Pkg.add("Ipopt")
Pkg.add("SpecialFunctions")
Pkg.add("ForwardDiff")
Pkg.add("AmplNLWriter")
Pkg.add("Cbc")
Pkg.add("JLD")


Pkg.add("IterTools")
Pkg.add("Tracker")
Pkg.add("Zygote")
Pkg.add("Surrogates")
Pkg.add("QuadGK")


using DataFrames, XLSX
using Missings
using Plots
using SpecialFunctions
using LinearAlgebra, Random, Distributions, Complementarity
using Test
using JuMP
using Ipopt

## load data
xf = XLSX.readxlsx("node_stats_forsimulation_all.xlsx") #xf["BHCs"]["A1:I1"]
data = vcat( [(XLSX.eachtablerow(xf[s]) |> DataFrames.DataFrame) for s in XLSX.sheetnames(xf)]... )
unique!(data) # delete duplicate rows
#data[(isequal.(data.tkr,"JPM") .| isequal.(data.tkr,"AIG")) .& isequal.(data.qt_dt,192),:]


#= network primitives from the data
p_bar = total liabilities
c = outside assets
assets = total assets
w = net worth
b = outside liabilities
=#

data = data[isless.(0,data.w), :]
data = data[isequal.(data.qt_dt,192), :] # keep quarter == 195 = 2008q4
sort!(data, :assets, rev = true)
data = data[1:5,:] # keep small number of nodes, for testing
N = size(data,1) # number of nodes

# rescale units
units = 1e6;
data[:,[:w, :c, :assets, :p_bar, :b]] .= data[!,[:w, :c, :assets, :p_bar, :b]]./units

# fake missing data to make problem feasible
data.b[:] .= missing
data.c[:] .= missing

# create network variables
#data.f = data.p_bar .- data.b # inside liabilities
#data.d =  data.assets .- data.c # inside assets

# keep track of missing variables
col_with_miss = names(data)[[any(ismissing.(col)) for col = eachcol(data)]] # columns with at least one missing 
#data_nm = dropmissing(data, disallowmissing=true) # drop missing
data_nm = coalesce.(data, 0.01) # replace missing by a value
nm_c = findall(x->x==0,ismissing.(data.c))
nm_b = findall(x->x==0,ismissing.(data.b))

# take a look
names(data) # column names
describe(data)

show(data, true)

# parameters
g0 = 0.0 # bankruptcy cost

# function shock(c)
#     # draw shocks that produce the probability of default δ 
#     a = 1
#     b = log.(data.delta)./(log.(1 .-data.w./c))
#     dist = Beta.(a,[b...])
#     draws = rand.(dist, 1)
#     vcat(draws...).*c
#     #rand(Float64,N) # normal
# end

## Optimization

# initial values

#rng = MersenneTwister(1234);
# dist = Beta.(2,fill(2,N,D))
# x0 = rand.(dist, 1)

# x0 = fill(0.0,N,D)
# A0 = rand(rng,N,N);[A0[i,i]=0.0 for i=1:N];A0=LowerTriangular(A0);


temp = convert(Array{Float32},data[:,[:delta, :delta_alt, :w, :assets, :p_bar]])
delta = copy(temp[:,1]); 
delta_alt = copy(temp[:,2]);
w= copy(temp[:,3]); 
assets= copy(temp[:,4]);
p_bar = copy(temp[:,5]);

g0 = 0.0f0   # bankruptcy cost
rng =  Random.seed!(123)
α=1.0f0 ;β=2.0f0 ;
dist = Beta(α,β)
dist = convert(Beta{Float32},dist)  
betapdf(x) = x.^(α .- 1) .* (1 .- x).^(β - 1) 
D = 1 # number of draws
x0 =  [0.4786082776942293
 0.6261009289132826
 0.0643506085692449
 0.32494449525126634
 0.5115660459579572];
# σA = convert(Array{Float32},rand(rng,N,N)) 
# σA[diagind(σA)] .= -1 
# for i=1:N
#     for j=i+1:N
#         σA[i,j] = -1 
#     end
# end

# using Flux.NNlib
# A0 = (NNlib.hardtanh.(σA).+1)./2
A0 = rand(N,N)

# set up optimization
# m = Model(with_optimizer(Ipopt.Optimizer, start_with_resto="yes", linear_solver="mumps",max_iter=10000,hessian_approximation="limited-memory"))
m = Model(with_optimizer(Ipopt.Optimizer))

@variable(m, 0<=p[i=1:N,j=1:D]<=data.p_bar[i], start = data.p_bar[i]) 
@variable(m, 0<=A[i=1:N, j=1:N]<=1, start=A0[i,j])  
@variable(m, 0<=c[i=1:N]<=data.assets[i], start = data_nm.c[i])  
@variable(m, 0<=b[i=1:N]<=data.p_bar[i], start = data_nm.b[i])   

@constraint(m, sum(A,dims=2).*data.p_bar .== data.p_bar .- b ) # payments to other nodes add up to inside liabilities f
@constraint(m, A' * data.p_bar .== data.assets .- c ) # payments from other nodes add up to inside assets d

# liabilities are net liabilities: A[i,i]=0 and A[i,j]A[j,i]=0
@constraint(m, [i = 1:N], A[i,i]==0)
for i=1:N
    j=1
    while j < i
        @complements(m, 0 <= A[i,j],  A[j,i] >= 0)
        j += 1
    end
end

# register max and min non-linear functions
maxfun(n1, n2) = max(n1, n2)
minfun(n1, n2) = min(n1, n2)

JuMP.register(m, :maxfun, 2, maxfun, autodiff=true)
JuMP.register(m, :minfun, 2, minfun, autodiff=true)

# clearing vector
myexpr = []
for j=1:D
    push!(myexpr, (1+g0).*(A'*p[:,j] .+ c .- x0[:,j].*c ) .- g0.*data.p_bar)
end
    
@variable(m, aux[i=1:N,j=1:D])
for i = 1:N
    for j=1:D
        @constraint(m, aux[i,j] .== myexpr[j][i] )
    end
end

for i = 1:N
    for j=1:D
        @NLconstraint(m,  minfun(data.p_bar[i], maxfun(aux[i,j],0)) == p[i,j] ) 
    end
end

# [fix(c[i], data.c[i]; force=true) for i  in nm_c]
# [fix(b[i], data.b[i]; force=true) for i  in nm_b]

#@NLobjective(m, Max , sum(sum( x0[i,j]*c[i]+data.p_bar[i]-p[i,j] for i=1:N)/D for j=1:D) ) #*sum(x[i]+p_bar[i]-p[i] for i=1:N) 
@NLobjective(m, Max , sum(sum( -p[i,j] for i=1:N)/D for j=1:D) ) 

#unset_silent(m)
JuMP.optimize!(m)

termination_status(m)
objective_value(m)

psol = JuMP.value.(p)
csol = JuMP.value.(c)
bsol = JuMP.value.(b)
Asol = JuMP.value.(A)
zsol = [Asol...;bsol...;csol...;psol...]
Zsol = zsol*zsol'
rank(Zsol)
eigvals(Zsol)

tol = 1e-6
@testset "check solution" begin
    @test norm( sum(Asol,dims=2).* data.p_bar .- (data.p_bar .- bsol)) < tol
    @test norm( Asol' * data.p_bar .- (data.assets .- csol)) < tol
   @test norm(diag(Asol)) < tol
   @test norm([Asol[i,j]*Asol[j,i] for i=1:N , j=1:N]) < tol
   @test all(0 .<=Asol.<=1)
    @test all(0 .<=bsol.<=data.p_bar)
    @test all(0 .<=csol.<=data.assets)
end

Aplot = copy(Asol)
Aplot[Aplot.<1f-4].=0

Sys.total_memory()/2^20 # total memory
Sys.free_memory()/2^20 # free memory
