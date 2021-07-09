import Pkg
Pkg.activate("benchmark_timing")
Pkg.instantiate()

using XLSX, DataFrames, Random, Test, NLsolve, Distributions, BenchmarkTools, LinearAlgebra, CUDA
using Surrogates, Flux
using Flux.Data: DataLoader
using IterTools: ncycle
using Flux: throttle
using ForwardDiff
using Convex, SCS, MathOptInterface, Parameters, Suppressor, Revise
using Dates, MonthlyDates
const MOI = MathOptInterface
using BSON
using BSON: @load, @save
include("NetworkUtils.jl")
using .NetworkUtils

# User Inputs
T=Float32
#N = 5 #number of ndoes
N = parse(Int64, ARGS[1])  #number of nodes from bach
Q = 200000 #number of samples - size(train) = Q/2 size(test) = Q/2
max_iter = 2*Q
dates=[Date(QuarterlyDate("2008-Q4"))]

# get network
# net = Network{T,N}()
data_dict = BSON.load("data.bson")
#net = netEmp(data_dict[:data],N)[1]
data_dict[:data].julia_dates = toJuliaDates(data_dict[:data].qt_dt)

keep = data_dict[:data].julia_dates.==lastdayofquarter(dates[1])
df_t = sort(data_dict[:data][keep,:], :assets, rev = true)
  
nm_b = [findall(x->x==false,ismissing.(df_t.b))...]
nm_b = nm_b[isless.(nm_b, N+1)] # indices for non-missing b
nm_c = [findall(x->x==false,ismissing.(df_t.c))...]
nm_c = nm_c[isless.(nm_c, N+1)] # indices for non-missing b
df_t = coalesce.(df_t, df_t.assets/1.5) # replace missing with some value
      
b = [deepcopy(Array{T}(df_t[1:N,[:b]]))...]
c = [deepcopy(Array{T}(df_t[1:N,[:c]]))...]
w = [deepcopy(Array{T}(df_t[1:N,[:w]]))...]
p_bar = [deepcopy(Array{T}(df_t[1:N,[:p_bar]]))...]
delta = [deepcopy(Array{T}(df_t[1:N,[:delta]]))...]      

# initial guess
rng =  Random.seed!(123)
A0 = zeros(N,N)
c0 =  c
b0 = b
α0 = fill(1.5,N)
β0= fill(15.0,N)

# Generate Inputs

# making X
function gen_x(N,Q,α0,β0)
    x = zeros(N,Q)
   for i = 1:N
        for j = 1:Q
            x[i,j] = rand(Beta(α0[i]+rand(1)[1]-0.5,β0[i]+ 12* (rand(1)[1]-0.5)), 1)[1] #drawing from distributions close to Beta(α0,β0)  but not quite. 
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
        A=rand(N,N)./(30 .* rand(1))
        R = 0.01* (s*rand(N, N).-s/2)
        for i = 1:N
            A[i,i] = 0.0
        end
        A_in = hcat(A_in, vec(A + R))
    end
    return A_in
end

halfQ =Q÷2
A_in_orig = mod(Q,2)==0 ? gen_A(N,halfQ) : gen_A(N,halfQ+1)
A_in_alt = NetworkUtils.trainingSetA(N,sample_size=halfQ,max_iter=max_iter,con="netting") # N=5 nodes, can use for training NN
A_in = hcat(A_in_orig,hcat(reshape.(A_in_alt,N^2)...))
# for i = 1:N^2
#     for j = 1:Int(Q/2)
#         A_in[i,j] = vcat(A_in_alt[j]...)[i]
#     end
# end

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
m = Int(Q÷2) #number of observations in test set
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

loss(x,y) = mean((model(x).-y).^2)
model = gpu(Chain(Dense(N*N+3*N, N*3, σ), Dense(N*3,N*2,σ), Dense(N*2,N*2,σ),Dense(N*2,N*2,σ), Dense(N*2,N,σ),Dense(N,N,σ), Dense(N, N)))
#@load "clearing_p_NN_gpu_N$N.bson" model
x_train = gpu(x_train)
y_train = gpu(y_train) 
x_test = gpu(x_test)
y_test = gpu(y_test) 

ps = Flux.params(model)
train_loader = gpu(DataLoader((x_train, y_train), batchsize=200, shuffle=true))

# Choosing Gradient Descent option
opt = ADAM(0.001, (0.9, 0.8))

## Setting up Call backs
evalcb() = @show(loss(x_test,y_test))
init_loss = loss(x_test,y_test)

## Training
Flux.@epochs 1500 Flux.train!(loss, ps, ncycle(train_loader, 5), opt, cb  = throttle(evalcb,150))

model = cpu(model)
@save "clearing_p_NN_gpu_N$N.bson" model
model = gpu(model)
## Testing
@testset "Neural Net Check" begin
    #Training set MSE and Abs Error
    @test Flux.mse(model(x_train),y_train) < 0.01 
    @test mean(abs, model(x_train) - y_train) < 0.01
    #Test set MSE and Abs Error
    @test Flux.mse(model(x_test),y_test) < 0.01  
    @test mean(abs, model(x_test) - y_test) < 0.01
    #Test set each prediction within 0.01
    pct_failing_train = sum(abs.(model(x_train) - y_train) .> 0.01)/(size(y_train,1)*size(y_train,2))
    pct_failing_test = sum(abs.(model(x_test) - y_test) .> 0.01)/(size(y_test,1)*size(y_test,2))
    @test pct_failing_train == 0
    @test pct_failing_test == 0
end;

#=

model = cpu(model) # move back to cpu to save
@save "clearing_p_NN_gpu.bson" model
#@load "clearing_p_NN_gpu.bson" model

model = gpu(model)
opt = ADAM(1e-6, (0.9, 0.8))
Flux.@epochs 1500 Flux.train!(loss, ps, ncycle(train_loader, 10), opt, cb  = throttle(evalcb,150))
model = cpu(model) # move back to cpu to save
@save "clearing_p_NN_gpu.bson" model

model = gpu(model)
opt = NADAM()
Flux.@epochs 150 Flux.train!(loss, ps, ncycle(train_loader, 10), opt, cb  = throttle(evalcb,150))
model = cpu(model) # move back to cpu to save
@save "clearing_p_NN_gpu.bson" model

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

=#
