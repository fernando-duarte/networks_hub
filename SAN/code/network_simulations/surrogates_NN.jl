import Pkg
Pkg.activate("benchmark_timing")
Pkg.instantiate()
Pkg.add("CUDA")
using XLSX, DataFrames, Random, Test, NLsolve, Distributions, BenchmarkTools, LinearAlgebra, CUDA
using Surrogates, Flux

# User Input
N = 15 # number of nodes
Q = 1  #number of beta draws

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
A=rand(N,N)./50
for i = 1:N
    A[i,i] = 0.0
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

x = x[:,1]
n = 1000
bounds = zeros(N), ones(N)
x_train = Surrogates.sample(n, bounds..., SobolSample())
y_train = zeros(N,n)
for i = 1:n
    y_train[:,i] = p_func(x_train[i],p_bar,A,c)
end

model = Chain(Dense(N, 20, relu), Dense(20,20,relu), Dense(20, N))
loss(x, y) = Flux.mse(model(x), y)

learning_rate = 0.1
optimizer = Descent(learning_rate)  # Simple gradient descent. See Flux documentation for other options.
n_epochs = 100
sgt = NeuralSurrogate(x_train, y_train, bounds..., model=model, loss=loss, opt=optimizer, n_echos=n_epochs)

x_test = Surrogates.sample(30, bounds..., SobolSample())
test_error = mean(abs2, sum(abs2, sgt(x) .- p_func(x,p_bar,A,c)) for x in x_test)

@btime sgt(x_test[1])
@btime p_func(x_test[1],p_bar,A,c)

## Sgt is faster than p_func, but face up front cost of training
