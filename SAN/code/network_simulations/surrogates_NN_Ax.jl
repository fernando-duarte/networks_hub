import Pkg
Pkg.activate("benchmark_timing")
Pkg.instantiate()
Pkg.add("CUDA")
Pkg.add("Flux")
Pkg.add("IterTools")
using XLSX, DataFrames, Random, Test, NLsolve, Distributions, BenchmarkTools, LinearAlgebra, CUDA
using Surrogates, Flux
using Flux.Data: DataLoader
using IterTools: ncycle

# User Inputs
N = 10 #number of points
Q = 1000 #number of samples

# Loading Data
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

x = zeros(N,Q)
for i = 1:N
    x[i,:] = rand(Beta(α0[i],β0[i]), Q)
end

# making concatination of many A matricies
A_in = rand(N,N)./50
for i = 1:N
    A_in[i,i] = 0.0
end
A_in = vec(A_in)
for i = 1:Q-1
    A=rand(N,N)./50
    for i = 1:N
        A[i,i] = 0.0
    end
    A_in = hcat(A_in, vec(A))
end
c = c0


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

x_train = vcat(x,A_in)
x_test = x_train[:,Q-29:end]
x_train = x_train[:,1:Q - 30]
y_train = zeros(N,Q-30)
y_test = zeros(N,30)
for i = 1:Q-30
    y_train[:,i] = p_func(x_train[1:N,i],p_bar,reshape(x_train[N+1:end,i], (N,N)),c) 
end
for i = 1:30
    y_test[:,i] = p_func(x_test[1:N,i],p_bar,reshape(x_test[N+1:end,i], (N,N)),c) 
end

model = Chain(Dense(N*N+N, 20, relu), Dense(20,20,relu), Dense(20, N))
loss(x, y) = Flux.mse(model(x), y)

ps = Flux.params(model)
train_loader = DataLoader((x_train, y_train), batchsize=2, shuffle=true)
opt = Descent()
Flux.@epochs 5 Flux.train!(loss, ps, ncycle(train_loader, 10), opt)

test_error = mean(sum(abs2, Array(model(x_test[:,i])) .- p_func(x_test[1:N,i],p_bar,reshape(x_train[N+1:end,i], (N,N)),c)) for i in 1:M)
