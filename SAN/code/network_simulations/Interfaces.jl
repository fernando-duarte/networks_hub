import Pkg
Pkg.activate("optimization")
Pkg.instantiate()

using Parameters, SparseArrays, Test
using Flux, BSON, DataFrames
using BSON: @load, @save
include("NetworkUtils.jl")
const NU = NetworkUtils
# using .NetworkUtils


## instantiate
dense = true
T=Float64
N=5
## network from DataFrame
data_dict = BSON.load("data.bson")
net = NU.netEmp(data_dict[:data],5)[1]

## optimization

# set up optimization
mini_batch = 1;
if dense
    param = NU.optimParam(net,mini_batch; dense=true)
else
    param = NU.optimParam(net,mini_batch)
end
@unpack M, vars = param

# z = vcat(A...,b,c,α,β,p...)

# upper and lower bounds
# low<=z<=upp
lowα = 0.5 ;uppα = 10
lowβ = 1   ;uppβ = 150
low_zero = eps(T) # eps(T) T(0.005)
@unpack N, p_bar, a = net
low = vcat(spzeros(T,N^2),fill(low_zero,N),fill(low_zero,N),fill(lowα,N),fill(lowβ,N),fill(low_zero,N*mini_batch))
upp = vcat(ones(T,N^2),p_bar,a,Array{T}(fill(uppα,N)),Array{T}(fill(uppβ,N)),repeat(p_bar,mini_batch))

# non-linear constraints cL<=c(z)<=cU, set cL=cU for equality constraints
cL = spzeros(3N+1+mini_batch*N)
cU = cL

# initial guess
if dense
    low = Array(low)
    cL = Array(cL)
    cU = Array(cU)
end
z0 = low + (upp-low).*rand(M)
vars0 = vars(z0)

##### ModelingToolkit #####
using ModelingToolkit, GalacticOptim, Flux, Optim, NonlinearSolve, LinearAlgebra


@variables z[1:2]
@parameters q[1:2]
eqs = [ 0 ~ z[1]*q[1] ,  0 ~ z[1]*q[1] + z[2]*q[2] + 3.4  ]
ns = NonlinearSystem(eqs, z, q)

guess = [ z[i] => Float64(i) for i=1:2]
ps =    [ q[1] => 1.0,  q[2] => 3.0]

prob = NonlinearProblem(ns,guess,ps; check_length=false)
sol = solve(prob,NewtonRaphson())

prob = NonlinearProblem(structural_simplify(ns),guess,ps)
sol = solve(prob,NewtonRaphson())

##### GalacticOptim #####

using GalacticOptim, Optim, Test, Random

rosenbrock(x, p) =  (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
x0 = zeros(2)
_p  = [1.0, 100.0]

function con2_c(x,p)
    [x[1]^2 + x[2]^2-1.0, x[2]*sin(x[1])-x[1]]
end

optprob = OptimizationFunction(rosenbrock, GalacticOptim.AutoForwardDiff();cons= con2_c)
prob = OptimizationProblem(optprob, x0, _p, lcons = [0,-Inf], ucons = [0,Inf])
sol = solve(prob, Optim.BFGS())

generate_gradient(optprob)

##### NLPModels #####
using ADNLPModels, NLPModels, NLPModelsIpopt

function consjac(nlp, x)
    return (cons(nlp, x),jac(nlp, x))
end

x0 = ones(N)/10
f(z) = NU.obj(z,x0,param)
f(z0)
#ADNLPModel(z->NU.c_lin(z,param), z0, low, upp, c, lcon, ucon)
nlp = ADNLPModel(z->0.0, z0, zeros(M),zeros(M), z->NU.c_NL(z,param),zeros(3*N+1),zeros(3*N+1))

f_k, nabf_k = objgrad(nlp, z0)
c_k, G_k = consjac(nlp, z0)
nab2f_k = hess(nlp, z0)

stats = ipopt(nlp)
