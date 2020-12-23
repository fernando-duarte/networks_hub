## load packages
import Pkg
Pkg.add("SpecialFunctions")
using SpecialFunctions


Pkg.add("DataFrames")
Pkg.add("XLSX")
Pkg.add("Missings")
Pkg.add("NLsolve")
Pkg.add("Complementarity")
Pkg.add("NLopt")
Pkg.add("JuMP")
Pkg.add("Ipopt")

Pkg.add("ForwardDiff")
#Pkg.build("SpecialFunctions")

Pkg.add("GraphRecipes")
Pkg.build("GR")

using GraphRecipes
using Plots

Pkg.add("Gadfly")
using Gadfly
spy(Asol)

Pkg.add("LightGraphs")
Pkg.add("SimpleWeightedGraphs")

using LightGraphs
using SimpleWeightedGraphs


using DataFrames, XLSX
using Missings
using ForwardDiff
using LinearAlgebra, Random, Distributions, NLsolve, Complementarity
using Test
using NLopt
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
data = data[1:10,:] # keep small number of nodes, for testing
N = size(data,1) # number of nodes

# rescale units
units = 1e6;
data[:,[:w, :c, :assets, :p_bar, :b]] .= data[!,[:w, :c, :assets, :p_bar, :b]]./units

# fake missing data to make problem feasible
data.b[:] .= missing
data.c[:] .= missing

# create network variables
data.f = data.p_bar .- data.b # inside liabilities
data.d =  data.assets .- data.c # inside assets

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

function shock(c)
    # draw shocks that produce the probability of default δ 
    a = 1
    b = log.(data.delta)./(log.(1 .-data.w./c))
    dist = Beta.(a,[b...])
    draws = rand.(dist, 1)
    vcat(draws...).*c
    #rand(Float64,N) # normal
end

## Optimization

# initial values
D = 2 # number of draws
rng = MersenneTwister(1234);
dist = Beta.(2,fill(2,N,D))
x0 = rand.(dist, 1)

x0 = fill(0.0,N,D)
A0 = rand(rng,N,N);[A0[i,i]=0.0 for i=1:N];A0=LowerTriangular(A0);

# set up optimization

#m = Model(Ipopt.Optimizer) # settings for the solver
m = Model(with_optimizer(Ipopt.Optimizer, start_with_resto="yes", linear_solver="mumps"))
#m = optimizer_with_attributes(Ipopt.Optimizer, "start_with_resto" => "yes", "linear_solver"=>"mumps")

@variable(m, 0<=p[i=1:N,j=1:D]<=data.p_bar[i], start = data.p_bar[i]) 
@variable(m, 0<=c[i=1:N]<=data.assets[i], start = data_nm.c[i])  
@variable(m, 0<=b[i=1:N]<=data.p_bar[i], start = data_nm.b[i])   
@variable(m, 0<=A[i=1:N, j=1:N]<=1, start=A0[i,j])  # start=A0[i,j]

@constraint(m, sum(A,dims=2).*data.p_bar .== data_nm.f) # payments to other nodes add up to inside liabilities f
@constraint(m, A' * data.p_bar .== data_nm.d) # payments from other nodes add up to inside assets d

#fix(v::VariableRef, value::Number; force::Bool = false)
#delete(model, con)

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

# chen-mangassarian smoothing function for max
# maxfun_cm(n1, n2,param=1000) = n1-n2 + log(1 + exp(-param*(n1-n2)))/param # converges to max when param \to \infty 
# minfun_cm(n1,n2,param=1/10000) = n1 - param* log(1+exp((n1-n2)/param)) # converges to min when param \to 0

# maxfun(n1,n2) = maxfun_cm(n1,n2,100) 
# minfun(n1,n2) = minfun_cm(n1,n2,1/1000)

# fisher-burmeister smoothing
# maxfun_fb(n1,n2,param=1/100000) = ( n1 + n2 + sqrt( (n1-n2)^2 + param^2 ) - param )/2  # converges to max as param \to 0

JuMP.register(m, :maxfun, 2, maxfun, autodiff=true)
JuMP.register(m, :minfun, 2, minfun, autodiff=true)

# clearing vector
myexpr = []
for j=1:D
    push!(myexpr, (1+g0)*(A'*p[:,j] .+ c .- x0[:,j].*c) .- g0.*data.p_bar)
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

[fix(c[i], data.c[i]; force=true) for i  in nm_c]
[fix(b[i], data.b[i]; force=true) for i  in nm_b]

@NLobjective(m, Max , sum(sum( x0[i,j]*c[i]+data.p_bar[i]-p[i,j] for i=1:N)/D for j=1:D) ) #*sum(x[i]+p_bar[i]-p[i] for i=1:N) 

#unset_silent(m)
JuMP.optimize!(m)

termination_status(m)
objective_value(m)

psol = JuMP.value.(p)
csol = JuMP.value.(c)
bsol = JuMP.value.(b)
Asol = JuMP.value.(A)

tol = 1e-6
if !(all(ismissing.(data.f)))
    @test norm( skipmissing(sum(Asol,dims=2).* data.p_bar .- data.f) ) < tol
end
if !(all(ismissing.(data.d)))
    @test norm( skipmissing(Asol' * data.p_bar .- data.d)) < tol
end
if !(all(ismissing.(data.c)))
    @test norm( skipmissing(psol .- min.(data.p_bar, max.((1+g0)*(Asol'*psol .+ data.c .- x0) .-g0.*data.p_bar,0))) ) < tol
end

@test norm(diag(Asol)) < tol
@test norm([Asol[i,j]*Asol[j,i] for i=1:N , j=1:N]) < tol
@test all(0 .<=psol)
@test all(0 .<=Asol.<=1)

Aplot = deepcopy(Asol)
Aplot[Aplot.<1e-3] .=0

# attributes here: https://docs.juliaplots.org/latest/generated/graph_attributes/
#method `:spectral`, `:sfdp`, `:circular`, `:shell`, `:stress`, `:spring`, `:tree`, `:buchheim`, `:arcdiagram` or `:chorddiagram`.

graphplot(LightGraphs.DiGraph(Aplot),
          nodeshape=:circle,
          markersize = 0.05,
          node_weights = data.assets,
          markercolor = range(colorant"yellow", stop=colorant"red", length=n),
          names = data.nm_short,
          fontsize = 8,
          linecolor = :darkgrey,
          edgewidth = (s,d,w)->500*Asol[s,d],
          arrow=true,
          method= :circular, #:chorddiagram,:circular,:shell
          )

# Pkg.add("JLD")
# using JLD
# save("/home/ec2-user/SageMaker/Test-AWS/net_opt.jld", "Asol", Asol,"data",data)



# using Flux, Zygote, Optim, FluxOptTools, Statistics
# m      = Chain(Dense(1,3,tanh) , Dense(3,1))
# x      = LinRange(-pi,pi,100)'
# y      = sin.(x)
# loss() = mean(abs2, m(x) .- y)
# Zygote.refresh()
# pars   = Flux.params(m)
# lossfun, gradfun, fg!, p0 = optfuns(loss, pars)
# res = Optim.optimize(Optim.only_fg!(fg!), p0, Optim.Options(iterations=1000, store_trace=true))


# import Pkg; Pkg.add("DistributedArrays")

# using Distributed
# using DistributedArrays
# #t = @async addprocs(8)
# t = addprocs(8)
# nprocs()
# nworkers()

# d=dzeros(100,100)
# d[:L]

# rmprocs(workers())

# @everywhere using SharedArrays


# # Utilize workers as and when they come online
# if nprocs() > 1   # Ensure at least one new worker is available
   
# end

# if istaskdone(t)   # Check if `addprocs` has completed to ensure `fetch` doesn't block
#     if nworkers() == N
#         new_pids = fetch(t)
#     else
#         fetch(t)
#     end
# end

# S = SharedArray{Int,2}((3,4), init = S -> S[localindices(S)] = repeat([myid()], length(localindices(S))))
# 3×4 SharedArray{Int64,2}:

# length(Sys.cpu_info())
# Threads.nthreads()


# addprocs(8)

# #BLAS.set_num_threads(1)
# @everywhere using <modulename> or @everywhere include("file.jl").


# Sys.free_memory()/2^20

#  versioninfo(verbose=true)



function contraction(p,x,c)
        min.(data_nm.p_bar, max.((1+g0)*(A0'*p .+ c .- x.*c) .- g0.*data_nm.p_bar,0)) 
end
contraction_iter(x, n::Integer,c) = n <= 0 ? data_nm.p_bar  : contraction(contraction_iter(x,n-1,c),x,c)

pinf = nlsolve(p->contraction(p, x0,data_nm.c)-p,p0, autodiff = :forward)
contraction_iter(x0,100,data_nm.c)-pinf.zero

using Zygote
Zygote.gradient(c -> -sum(contraction_iter(x0,4,c)), data_nm.c)

a = 1
b = log.(data_nm.delta)./(log.(1 .- data_nm.w./data_nm.c))
dist = Beta.(a,[b...])

histogram(pdf.(dist,rand(10)))

shock() = vcat(rand.(dist, 1)...)
x0 = shock()
Zygote.gradient(c -> -sum(contraction_iter(x0,4,c)), data_nm.c)
function loss(c)
    x0 = shock()
    sum(contraction_iter(x0,4,c))
end
Zygote.gradient(c -> loss(c), data_nm.c)
check DistributionsAD.jl, Turing.jl
Flux.train!
using CUDA

# Test GPU movement inside the call to `gradient`
@testset "GPU movement" begin
  r = rand(Float32, 3,3)
  @test gradient(x -> sum(cu(x)), r)[1] isa AbstractArray
end

@testset "basic bcasting" begin
  a = cu(Float32.(1:9))
  v(x, n) = x .^ n
  pow_grada = cu(Float32[7.0, 448.0, 5103.0, 28672.0, 109375.0, 326592.0, 823543.0, 1.835008e6, 3.720087e6])
  @test gradient(x -> v(x, 7) |> sum, a) == (pow_grada,)
  w(x) = broadcast(log, x)
  log_grada = cu(Float32[1.0, 0.5, 0.33333334, 0.25, 0.2, 0.16666667, 0.14285715, 0.125, 0.11111111])
  @test gradient(x -> w(x) |> sum, a) == (log_grada,)
end


# Zygote.hessian(c1 -> -sum(contraction_iter(x0,4,[c1 data_nm.c[2:end]])), data_nm.c)

# Zygote.hessian(((a, b),) -> a*b, [2, 3])

# using Zygote: @adjoint
# nestlevel() = 0
# @adjoint nestlevel() = nestlevel()+1, _ -> nothing
# function f(x)
#     println(nestlevel(), " levels of nesting")
#     return x^2
# end
# grad(f, x) = gradient(f, x)[1]
# grad(f,1)
