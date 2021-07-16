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
N = 5 #number of points
Q = 100000 #number of samples - size(train) = Q/2 size(test) = Q/2
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

function p_func_true(x)
    sum(p_func(x[N*N+2*N+1:end],p_bar,reshape(collect(x[1:N*N]), (N,N)),x[N*N+N+1:N*N+2*N]))
end

## Training NN
m = Int(Q÷2) #number of observations in test set
x_train = vcat(A_in, b_in, c_in, x_in) #order is A, b,c,x
x_train = x_train[:,randperm(Q)] #randomly sort columns of X_train
x_test = x_train[:,1:m]
x_train = x_train[:,m+1:end]
y_train = zeros(1,Q-m)
y_test = zeros(1,m)
for i = 1:Q-m
    y_train[i] = p_func_true(x_train[:,i])
end
for i = 1:m
    y_test[i] = p_func_true(x_test[:,i])
end

loss(x,y) = mean((model(x).-y).^2)
x_train = gpu(x_train)
y_train = gpu(y_train) 
x_test = gpu(x_test)
y_test = gpu(y_test) 

model = gpu(Chain(Dense(N*N+3*N, N^2, σ), Dense(N^2,N*3,σ), Dense(N*3,N*2,σ), Dense(N*2, N*2, σ), Dense(N*2,N*2,σ),Dense(N*2,N*2,σ), Dense(N*2,N,σ),Dense(N,N,σ), Dense(N, 1)))
ps = Flux.params(model)
train_loader = gpu(DataLoader((x_train, y_train), batchsize=500, shuffle=true))

## Setting up Call backs
evalcb() = @show(loss(x_train,y_train), loss(x_test,y_test))
init_loss = loss(x_test,y_test)

## Training
opt = ADAM(0.01, (0.9, 0.8))
Flux.@epochs 100 Flux.train!(loss, ps, ncycle(train_loader, 5), opt, cb  = throttle(evalcb,150))
opt = ADAM(0.001, (0.9, 0.8))
Flux.@epochs 500 Flux.train!(loss, ps, ncycle(train_loader, 5), opt, cb  = throttle(evalcb,150))
opt = ADAM(0.0001, (0.9, 0.8))
Flux.@epochs 100 Flux.train!(loss, ps, ncycle(train_loader, 5), opt, cb  = throttle(evalcb,150))
opt = ADAM(0.00005, (0.9, 0.8))
Flux.@epochs 100 Flux.train!(loss, ps, ncycle(train_loader, 5), opt, cb  = throttle(evalcb,150))

model = cpu(model)
@save "clearing_p_NN_gpu_sum_N$N.bson" model
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
