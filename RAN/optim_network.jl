## load packages
using Pkg
Pkg.add("DataFrames")
Pkg.add("XLSX")
Pkg.add("Missings")

using DataFrames, XLSX
using Missings

Pkg.add("NLsolve")
Pkg.add("JuMP")
Pkg.add("Ipopt")
Pkg.add("Complementarity")
Pkg.add("NLopt")

using LinearAlgebra, Random, Distributions, NLsolve
using Test
using NLopt
using JuMP, Ipopt
using Complementarity


## load data
xf = XLSX.readxlsx("node_stats_for_simulation.xlsx") #xf["BHCs"]["A1:I1"]
data = vcat( [(XLSX.eachtablerow(xf[s]) |> DataFrames.DataFrame) for s in XLSX.sheetnames(xf)]... )
unique!(data) # delete duplicate rows
#sort!(data, :assets, rev = true)
data = data[1:30,:] # keep small number of nodes, for testing

#= network primitives from the data
p_bar = total liabilities
c = outside assets
assets = total assets
w = net worth
b = outside liabilities
=#

data = data[data.w .> 0, :]
N = size(data,1) # number of nodes

# rescale units
units = 1e6;
data[:,[:w, :c, :assets, :p_bar, :b]] .= data[!,[:w, :c, :assets, :p_bar, :b]]./units

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
x0 = 0.0*similar(data.p_bar)

rng = MersenneTwister(1234);
A0 = rand(rng,N,N);[A0[i,i]=0.0 for i=1:N];A0=LowerTriangular(A0);

# set up optimization

m = Model(Ipopt.Optimizer) # settings for the solver

@variable(m, 0<=p[i=1:N]<=data.p_bar[i], start = data.p_bar[i])
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
myexpr = (1+g0)*(A'*p .+ c .- x0) .- g0.*data.p_bar
@variable(m, aux[i=1:N])
@constraint(m, aux .== myexpr )
for i = 1:N
    @NLconstraint(m,  minfun(data.p_bar[i], maxfun(aux[i],0)) == p[i] )
    #@NLconstraint(m,  min(p_bar0[i], max(aux[i],0)) == p[i] )
end

[fix(c[i], data.c[i]; force=true) for i  in nm_c]
[fix(b[i], data.b[i]; force=true) for i  in nm_b]

@NLobjective(m, Max , sum( x0[i]+data.p_bar[i]-p[i] for i=1:N) ) #*sum(x[i]+p_bar[i]-p[i] for i=1:N)

unset_silent(m)
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


Pkg.add("JLD")
using JLD
save("/home/ec2-user/SageMaker/Test-AWS/net_opt.jld", "Asol", Asol,"data",data)



# using Flux, Zygote, Optim, FluxOptTools, Statistics
# m      = Chain(Dense(1,3,tanh) , Dense(3,1))
# x      = LinRange(-pi,pi,100)'
# y      = sin.(x)
# loss() = mean(abs2, m(x) .- y)
# Zygote.refresh()
# pars   = Flux.params(m)
# lossfun, gradfun, fg!, p0 = optfuns(loss, pars)
# res = Optim.optimize(Optim.only_fg!(fg!), p0, Optim.Options(iterations=1000, store_trace=true))
