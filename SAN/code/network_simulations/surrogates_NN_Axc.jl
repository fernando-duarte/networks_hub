import Pkg
Pkg.activate("benchmark_timing")
Pkg.instantiate()

using XLSX, DataFrames, Random, Test, NLsolve, Distributions, BenchmarkTools, LinearAlgebra, CUDA
using Surrogates, Flux
using Flux.Data: DataLoader
using IterTools: ncycle
using Flux: throttle
using ForwardDiff
using Convex, SCS, LinearAlgebra, Distributions, Test, Random, MathOptInterface, Parameters, Suppressor
const MOI = MathOptInterface
include("generateA.jl")
using .generateA


# User Inputs
N = 5 #number of points
Q = 1000 #number of samples
max_iter = 10000

# Loading Data
xf = XLSX.readxlsx("node_stats_forsimulation_all.xlsx") 
data = vcat( [(XLSX.eachtablerow(xf[s]) |> DataFrames.DataFrame) for s in XLSX.sheetnames(xf)]... ) #for s in XLSX.sheetnames(xf) if (s!="Aggregates Composition" && s!="Dealer Aggregates" && s!="Approx Aggregates")
unique!(data) # delete duplicate rows, use `nonunique(data)` to see if there are any duplicates
data = data[isequal.(data.qt_dt,195), :] # keep quarter == 195 = 2008q4
sort!(data, :assets, rev = true)
data = data[1:N,:] # keep small number of nodes, for testing
units = 1e6;
data[:,[:w, :c, :assets, :p_bar, :b]] .= data[!,[:w, :c, :assets, :p_bar, :b]]./units


col_with_miss = names(data)[[any(ismissing.(col)) for col = eachcol(data)]] # columns with at least one missing
data_nm = coalesce.(data, data.assets/1.5) # replace missing by a value
nm_c = findall(x->x==0,ismissing.(data.c))
nm_b = findall(x->x==0,ismissing.(data.b))
dropmissing(data, [:delta, :delta_alt, :w, :assets, :p_bar]) # remove type missing

names(data) # column names
describe(data)
show(data, allcols = true)

p_bar = data.p_bar
assets = data.assets
delta = data.delta
w= data.w
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
        val = sol.zero[1]
    else
        try
            sol = nlsolve(ff!,[2.0])
            val = sol.zero[1]
        catch
            val = 25.0
        end
    end
    push!(β0,val)
end

# Generate Inputs

# making X
function gen_x(N,Q,α0,β0)
    x = zeros(N,Q)
   for i = 1:N
        for j = 1:Q
            x[i,j] = rand(Beta(α0[i]+rand(1)[1]-0.5,β0[i]+ 10* (rand(1)[1]-0.5)), 1)[1] #drawing from distributions close to Beta(α0,β0)  but not quite. 
        end
    end
    return x
end
x_in = gen_x(N,Q,α0,β0)

# making A
function gen_A(N,Q)
    s =2 
    R = 0.01* (s*rand(N, N).-s/2)
    A_in = rand(N,N)./(50 .* rand(1))
    A_in = A_in + R
    for i = 1:N
        A_in[i,i] = 0.0
    end
    A_in = vec(A_in)
    for i = 1:Q-1
        A=rand(N,N)./(50 .* rand(1))
        R = 0.01* (s*rand(N, N).-s/2)
        for i = 1:N
            A[i,i] = 0.0
        end
        A_in = hcat(A_in, vec(A + R))
    end
    return A_in
end

A_in = gen_A(N,Q)
A_in_alt = trainingSetA(N,sample_size=Int(Q/2),max_iter=max_iter,con="netting") # N=5 nodes, can use for training NN
for i = 1:N^2
    for j = 1:Int(Q/2)
        A_in[i,j] = vcat(A_in_alt[j]...)[i]
    end
end

# Making c
function gen_c(N,Q)
    c_in = zeros(N,Q)
    for i = 1:Q
        c_in[:,i] = c0 .+ ((rand(N) .- 0.5)./5) #slight permutations to c
    end
    ind = c_in .< 0 
    c_in[ind] .= 0.0
    return c_in
end
c_in = gen_c(N,Q)


# Making b
function gen_b(N,Q)
    b_in = zeros(N,Q)
    for i = 1:Q
        b_in[:,i] = b0 .+ ((rand(N) .- 0.5)./5) #slight permutations to c
    end
    ind = b_in .< 0 
    b_in[ind] .= 0.0
    return b_in
end
b_in = gen_b(N,Q)

# Function to replicate
function clearing_vec!(F, p, x, p_bar, A, c)
    # when the system is solved, F is zero
    F .= p.-min.(p_bar,max.(A'*p .+c .- x,0))
end

# clearing vector as a function of the shock x
function p_func(x, p_bar, A, c)
  sol = nlsolve((F, p)->clearing_vec!(F, p, x, p_bar, A, c),[p_bar...], autodiff = :forward)
  return sol.zero; # the solution to the system
end

## Training NN
m = Int(Q/10) #number of observations in test set
x_train = vcat(A_in, b_in, c_in, x_in) #order is A, b,c,x
x_train = x_train[:,randperm(Q)] #randomly sort columns of X_train
x_test = x_train[:,1:m]
x_train = x_train[:,m+1:end]
y_train = zeros(N,Q-m)
y_test = zeros(N,m)
for i = 1:Q-m
    y_train[:,i] = p_func(x_train[N*N+2*N+1:end,i],p_bar,reshape(x_train[1:N*N,i], (N,N)),x_train[N*N+N+1:N*N+2*N,i])  #x, p_bar, A, c
end
for i = 1:m
    y_test[:,i] = p_func(x_test[N*N+2*N+1:end,i],p_bar,reshape(x_test[1:N*N,i], (N,N)),x_test[N*N+N+1:N*N+2*N,i]) 
end

model = Chain(Dense(N*N+3*N, 50, σ), Dense(50,20,σ), Dense(20,10,σ), Dense(10,5,σ), Dense(5, N))
loss(x, y) = Flux.mse(model(x), y)

ps = Flux.params(model) #parameters that update
int_params = ps
train_loader = DataLoader((x_train, y_train), batchsize=100, shuffle=true)

# Choosing Gradient Descent option
opt = ADAM(0.001, (0.9, 0.8))

## Setting up Call backs
evalcb() = @show(loss(x_test,y_test))
init_loss = loss(x_test,y_test)

## Training
Flux.@epochs 150 Flux.train!(loss, ps, ncycle(train_loader, 10), opt, cb  = throttle(evalcb,50))

# Computing Jacobian and Hessian
j = x -> ForwardDiff.jacobian(model, x)
j(x_test[:,1])

## Creating many functions for each output of the NN
function model_hessian(x) 
    sum(model(x))
end

h_model = x -> ForwardDiff.hessian(model_hessian, x)
h_model(x_test[:,1])

function p_func2(w)
    x = w[N*N+2*N+1:end]
    A = reshape(w[1:N*N],(N,N))
    c = w[N*N+N+1:N*N+2*N]
  sol = nlsolve((F, p)->clearing_vec!(F, p, x, p_bar, A, c),[p_bar...], autodiff = :forward)
  return sol.zero; # the solution to the system
end

function p_hessian(x) 
    sum(p_func2(x))
end

h_p = x -> ForwardDiff.hessian(p_hessian, x)
p_hessian(x_test[:,1])
h_p(x_test[:,1]) #doesn't work unclear why

@testset "Neural Net Check" begin
    W = 100 #size of test set
    c_in = gen_c(N,W)
    b_in = gen_b(N,W)
    x_in = gen_x(N,W,α0,β0)
    
    A_in = gen_A(N,W)
    A_in_alt = trainingSetA(N,sample_size=Int(W/2),max_iter=max_iter,con="netting") # N=5 nodes, can use for training NN
    for i = 1:N^2
        for j = 1:Int(W/2)
           A_in[i,j] = vcat(A_in_alt[j]...)[i]
        end
    end
    x_test = vcat(A_in, b_in, c_in, x_in) #order is A, b,c,x
    x_test = x_test[:,randperm(W)] #randomly sort columns of X_train
    y_test = zeros(N,W)
    for i = 1:W
        y_test[:,i] = p_func(x_test[N*N+2*N+1:end,i],p_bar,reshape(x_test[1:N*N,i], (N,N)),x_test[N*N+N+1:N*N+2*N,i]) 
    end

    @test loss(x_test,y_test) < 0.01
      
end;
