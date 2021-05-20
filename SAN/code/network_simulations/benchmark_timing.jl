# Found I have to do pkg.add() for the first run
import Pkg
Pkg.activate("joint_timing")
Pkg.instantiate()

using Test, BenchmarkTools, Profile
using Random, SpecialFunctions, LinearAlgebra, ForwardDiff, FiniteDiff, Zygote, Distributions, DistributionsAD
using JuMP, Ipopt, NLsolve
using DataFrames, XLSX, Missings, JLD2, CSV
using Quadrature, Cuba, Cubature, HCubature

display("Finished setting up packages")

# User Inputs
profile = 1 # indicator if you want to profile the function  (1 if so)
max_iters = 10000 #maximum number of iterations when optimizing
#N = parse(Int64, ARGS[1]) # number of nodesif called via bash 
N = 5

## load data
xf = XLSX.readxlsx("node_stats_forsimulation_all.xlsx") 
data = vcat( [(XLSX.eachtablerow(xf[s]) |> DataFrames.DataFrame) for s in XLSX.sheetnames(xf)]... ) #for s in XLSX.sheetnames(xf) if (s!="Aggregates Composition" && s!="Dealer Aggregates" && s!="Approx Aggregates")
unique!(data) # delete duplicate rows, use `nonunique(data)` to see if there are any duplicates
data = data[isequal.(data.qt_dt,195), :] # keep quarter == 195 = 2008q4
sort!(data, :assets, rev = true)
data = data[1:N,:] # keep small number of nodes, for testing
#N = size(data,1) # number of nodes
units = 1e6;
data[:,[:w, :c, :assets, :p_bar, :b]] .= data[!,[:w, :c, :assets, :p_bar, :b]]./units
# data.b[:] .= missing
# data.c[:] .= missing

col_with_miss = names(data)[[any(ismissing.(col)) for col = eachcol(data)]] # columns with at least one missing
data_nm = coalesce.(data, data.assets/1.5) # replace missing by a value
nm_c = findall(x->x==0,ismissing.(data.c))
nm_b = findall(x->x==0,ismissing.(data.b))
dropmissing(data, [:delta, :delta_alt, :w, :assets, :p_bar]) # remove type missing

names(data) # column names
describe(data)
show(data, allcols = true)

p_bar = convert(Array{Float64},data.p_bar)
assets = convert(Array{Float64},data.assets)
delta = convert(Array{Float64},data.delta)
w=convert(Array{Float64},data.w)
g0 = 0.05   # bankruptcy cost

# initial guess
rng =  Random.seed!(123)
A0 = zeros(N,N)
c0 =  data_nm.c 
b0 = data_nm.b
α0 = 1.0*ones(N)
β0= []
for i=1:N
    function ff!(F, x)
        F[1] = 1-cdf(Beta(α0[i],x[1]),w[i]/c0[i])-delta[i] 
    end
    if i==1
        sol = nlsolve(ff!,[50.0])
    else
        sol = nlsolve(ff!,[1.0])
    end
    @test converged(sol)
    push!(β0,sol.zero[1])
end

dist_pdf(x,p) = exp(DistributionsAD.logpdf(DistributionsAD.Beta(p[1],p[2]),x[1]))

function dist_cdf(p...)
    α,β,wc = p[1], p[2], p[3]
    s(x,wc) = 1 ./(1 .+exp.(-1000*(x-wc))) # approximates indicator function of x>=wc
    obj(x,p) = dist_pdf(x,p)*s(x,wc...)
    prob = QuadratureProblem(obj,0.0,1.0,[α...,β...])
    Quadrature.solve(prob,QuadGKJL(),reltol=1e-5,abstol=1e-5)[1]
end

function dist_pdf_mv(x,params)
 # α,β,A,b,c = params[1][1:N], params[1][N+1:2*N], params[1][2*N+1:2*N+N^2], params[1][2*N+N^2+1:3*N+N^2], params[1][3*N+N^2+1:4*N+N^2]
     α,β,A,b,c = params[1:N], params[N+1:2*N], params[2*N+1:2*N+N^2], params[2*N+N^2+1:3*N+N^2], params[3*N+N^2+1:4*N+N^2]
     A = reshape([A...],N,N) 
     joint_pdf = 1.0;
     for i=1:N
         joint_pdf = joint_pdf*dist_pdf(x[i],[α[i],β[i]]) 
     end
 
    function contraction(x,p)
        min.(p_bar, max.((1 .+g0).*(A'*p .+ c .- x.*c) .- g0.*p_bar,0))
    end
    contraction_iter(x, n::Integer) = n <= 0 ? p_bar  : contraction(x,contraction_iter(x,n-1))
    loss_x(x) = -sum(contraction_iter(x,2)).*dist_pdf(x,[α...,β...])   
    loss_x(x)*joint_pdf
end

free_mem = []

function ev(params...)
    prob = QuadratureProblem(dist_pdf_mv,zeros(N),ones(N),[params...])
    global free_mem = [free_mem; Sys.free_memory()/1000/2^20]
    Quadrature.solve(prob,CubaSUAVE(),reltol=1e-5,abstol=1e-5)[1]
end

# set up optimization
m = Model(optimizer_with_attributes(Ipopt.Optimizer,"max_iter"=>max_iters,"print_level"=>5))

@variable(m, 0.8<=α[i=1:N]<=3.0, start = α0[i]) 
@variable(m, 1.0<=β[i=1:N]<=150.0, start = β0[i]) 
#@variable(m, 0<=p[i=1:N,j=1:D]<=p_bar[i], start = p_bar[i]) 
@variable(m, 0<=A[i=1:N, j=1:N]<=1, start=A0[i,j])  
@variable(m, w[i]<=c[i=1:N]<=assets[i], start = c0[i])  
@variable(m, 0<=b[i=1:N]<=p_bar[i], start = b0[i])   

@constraint(m, -sum(A,dims=2).*p_bar .+ p_bar .==  b ) # payments to other nodes add up to inside liabilities f
@constraint(m, A' * p_bar .== assets .- c ) # payments from other nodes add up to inside assets d
# liabilities are net liabilities: A[i,i]=0 and A[i,j]A[j,i]=0
@constraint(m, [i = 1:N], A[i,i]==0)
for i=1:N
    j=1
    while j < i
       @constraint(m, 0 == A[i,j] * A[j,i])
        j += 1
    end
end

# register max and min non-linear functions
maxfun(n1, n2) = max(n1, n2)
minfun(n1, n2) = min(n1, n2)
JuMP.register(m, :maxfun, 2, maxfun, autodiff=true)
JuMP.register(m, :minfun, 2, minfun, autodiff=true)

# match probabilities of default
JuMP.register(m, :dist_cdf, 3, dist_cdf, autodiff=true)
@NLexpression(m,wc[i=1:N],w[i]/c[i])
for i in 1:N
    @eval ($(Symbol("vv$i"))) = [α[$i],β[$i],wc[$i]]
    @eval @NLconstraint(m,dist_cdf(($(Symbol("vv$i")))...)== delta[$i]) 
end

[fix(c[i], data_nm.c[i]; force=true) for i in nm_c] #fixing c to data
#[fix(b[i], data.b[i]; force=true) for i  in nm_b] #fixing b to data
JuMP.register(m, :ev, 4*N+N^2, ev, autodiff=true)
vars = [α...,β...,A...,b...,c...]
@NLobjective(m, Min, ev(vars...)  )

unset_silent(m)


time_start = time()

#================= 
@benchmark and @btime run JuMP.optimize!(m) 4 or more times
@time and @elapsed run it once
Compare screen output when running:
```
@benchmark JuMP.optimize!(m) # 1 iter=  10.358 seconds, 3.95 total gib #as many iters as desired = 353.666 s, 133.93 gib
@btime JuMP.optimize!(m)  # 10.395 seconds, 3.95 total gib as many int as desired 347s, 133.93 gb
@time JuMP.optimize!(m)  # 10.342 seconds, 3.952 gib, 350 seconds, 133.927 gib
@elapsed JuMP.optimize!(m) # 10.38. 349.5 seconds
```
================#

# running optimization
if profile == 1
    Profile.clear() # unnecessary for single batch runs?
    Profile.init(delay=0.1) # delay per snapshot in secods, default is 0.001 -- increase for faster and coarser profiling
    mem_usage = @profile @allocated JuMP.optimize!(m)  
    time_usage = @profile @elapsed JuMP.optimize!(m)
else 
    mem_usage = @allocated JuMP.optimize!(m) 
    time_usage = @elapsed JuMP.optimize!(m)  
end

## recording meta_data
time_end = time()
total_mem = mem_usage/1000/2^20
total_time = time_end - time_start # ≈time_usage?

max_mem = max(free_mem...) - min(free_mem...)

st0 = termination_status(m)
obj0 = objective_value(m)

csol0 = JuMP.value.(c)
bsol0 = JuMP.value.(b)
Asol0 = JuMP.value.(A)
αsol0 = JuMP.value.(α)
βsol0 = JuMP.value.(β)
objective = objective_value(m)

cdf_pkg = [1-cdf(Beta(αsol0[i],βsol0[i]),w[i]/csol0[i]) for i=1:N]
cdf_opt = [dist_cdf([αsol0[i],βsol0[i],w[i]/csol0[i] ]...) for i=1:N]
    
## saving meta data

test_passed = (norm( sum(Asol0,dims=2).* data.p_bar .- (data.p_bar .- bsol0)) < 1e-5) + 
(norm( Asol0' * data.p_bar .- (data.assets .- csol0)) < 1e-5) + 
(norm(diag(Asol0)) < 1e-5) + (norm([Asol0[i,j]*Asol0[j,i] for i=1:N , j=1:N]) < 1e-5) + 
(all(0 .<=Asol0.<=1)) + (all(0 .<=bsol0.<=data.p_bar)) + (all(0 .<=csol0.<=data.assets)) +
(norm(cdf_pkg-cdf_opt)<0.01) + (norm(delta-cdf_opt)<0.01)

meta_data = DataFrame(nodes = N, total_time_min = total_time/60, total_mem_gb = total_mem, max_mem = max_mem, optimization_output = st0, tests_passed = test_passed, total_tests = 9, max_iters = max_iters)
CSV.write("meta_data_$(N)_SUAVE.csv", meta_data)

## Exporting Results to csv
data_out = DataFrame(node = 1:N, csol0 = csol0, bsol0 = bsol0, alphasol0 = αsol0, betasol0 = βsol0)
data_asol0 = convert(DataFrame, Asol0)
data_out = hcat(data_out, data_asol0)
data_out.objective = objective * ones(N)
CSV.write("data_out_$(N)_SUAVE.csv",  data_out)

# Saving Variables
@save "variables_$N.jld2" A0 Asol0 N b0 bsol0 c0 csol0 delta p_bar w α0 αsol0 β0 βsol0 assets g0 rng

# Saving profile output
if profile == 1
    open("profile_$N.txt", "w") do s
        Profile.print(IOContext(s, :displaysize => (24, 500)))
    end
end

if profile == 0
    open("profile_$N.txt", "w") do io
       println(io,"No Profiling Requested")
    end
end