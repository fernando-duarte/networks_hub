import Pkg
Pkg.activate("optimization")
Pkg.instantiate()

using Parameters, SparseArrays, Test
using Flux, BSON, DataFrames
using BSON: @load, @save
include("NetworkUtils.jl")
# using .NetworkUtils

# send negatives to zero
# send smallest of aij and aji to zero
# divide each element by row sum OR send aij>1 to 1, then divide all aij (OR just the ones moved to 1) by row sum (OR so that row sum =1 )

## instantiate
dense = true
T=Float64
N=5
## network from DataFrame
data_dict = BSON.load("data.bson")
net = netEmp(data_dict[:data],5)[1]

## optimization

# set up optimization
mini_batch = 1;
if dense
    param = optimParam(net,mini_batch; dense=true)
else
    param = optimParam(net,mini_batch)
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

# clearing vector using NN
#@load "clearing_p_NN_gpu.bson" model

using ForwardDiff, LinearAlgebra


x0 = rand(N)
p_guess = vars0.p
p_guess_rep = repeat(p_guess,1,mini_batch)
p_bar_rep = repeat(net.p_bar,1,mini_batch)
c_rep = repeat(vars0.c,1,mini_batch)
A = LowerTriangular(ones(N,N))/3
function clearing_p!(F, p, x, Aᵀ, c,γ,p_bar_rep,c_rep)
    F .= p.-min.(p_bar_rep,max.( (1+γ)*(Aᵀ*p .+ (one(x[1]) .- x).*c_rep) .- γ.*p_bar_rep, 0 ) )
end
F = ones(size(p_bar_rep))
p0 = deepcopy(p_bar_rep)
clearing_p!(F, p0, x0, Aᵀ, c,γ,p_bar_rep,c_rep)

s=nlsolve((F, p)->clearing_p!(F, p, x, transpose(A), vars0.c,net.γ,p_bar_rep,c_rep),p_bar_rep, autodiff = :forward)
pnl=s.zero
pnn=p_NN(z0,x0)
[p_bar,pnl,pnn]


using Nonconvex
# include("nlsolve_rules.jl")
# include("nlsolve_rrule.jl")
m_c = Nonconvex.Model(z->obj(z,x0,param; use_nlsolve=true))
Nonconvex.addvar!(m_c,low,upp)
Nonconvex.add_eq_constraint!(m_c,z->c_NL(z,param))
alg = IpoptAlg()
options = Nonconvex.IpoptOptions()
obj(z0,x0,param)
obj(z0,x0,param; use_nlsolve=true)

r_c = Nonconvex.optimize(m_c, alg,z0, options = options);


tol = 1e-12
@test sum(abs2,c_NL(r_c.minimizer,x0,param)) < tol



alg = MMA87() # or MMA02()
options = MMAOptions()



# empty constructor of Network defaults to an example from Glasserman and Young
# net = Network{T,5}()
# p_net = plot(net,plot_outside = true)

# # network from xlsx file
#fileName="node_stats_forsimulation_all.xlsx"
#net = netEmp(fileName,N)[1]
    
# fileName="glasserman_young_example.xlsx"
# N = 5 
# netGY2 = netEmp(fileName,N)[1]
# # TO DO: compare netGY with netGY2 

# plots
# pGY = plot(netGY)
# pGY_in = plot(netGY,plot_outside = false)
# pGY_out = plot(netGY,plot_outside = true)
# savefig(pGY, "GYnet.pdf") 
# savefig(pGY_in, "GYnet_in.pdf") 
# savefig(pGY_out, "GYnet_out.pdf") 

# pGY2 = plot(netGY2)
# savefig(pGY2, "GYnet.pdf") 

# pnet1=plot(net1)
# savefig(pnet1, "net.pdf") 

# ## create network
# T=Float64
# N=5
# net = Network{T,5}()

using NLsolve
using ForwardDiff, LinearAlgebra, Zygote

Aᵀ= (LowerTriangular(ones(N,N))-Matrix(I,N,N))/3
A = transpose(Aᵀ)
@unpack c,γ,p_bar = net
p_bar_rep = repeat(p_bar,1,mini_batch)
c_rep = repeat(c,1,mini_batch)
x = ones(N)/10
function contraction4!(F::AbstractArray{T}, p::AbstractArray{T}, x, Aᵀ, c,γ,p_bar_rep,c_rep) where T<:Real
    F .= min.(p_bar_rep,max.( (1+γ)*(Aᵀ*p + (one(x[1]) .- x).*c_rep) - γ.*p_bar_rep, 0 ) )
end
F0 = zeros(N,mini_batch)
p0 = deepcopy(p_bar_rep)
contraction4!(F0, p0, x, Aᵀ, c,γ,p_bar_rep,c_rep)
ForwardDiff.gradient(p->contraction4!(F0, p, x, Aᵀ, c,γ,p_bar_rep,c_rep),p0)


function contraction(p)
    min.(2,max.(p, 0 ) )
end
contraction(p0)
c2! = make_inplace(contraction)
ForwardDiff.jacobian(contraction,p0)
ForwardDiff.jacobian(c2!,p0)


ForwardDiff.derivative(y->max(x0,y),2)



# Function converter for out of place --> in place.
function make_inplace(f::Function)
    function inplace_f(out, x)
        out .= f(x)
    end
    return inplace_f
end


contraction_iter!(F, x, n::Integer) = n <= 0 ? F  : contraction!(F,contraction_iter!(F,x,n-1))
F= copy(p_bar_rep)
contraction_iter!(F,x0,10)
p_bar_rep

using ForwardDiff, Zygote
ForwardDiff.jacobian(p->contraction5!(F,p) ,p_bar_rep)


ForwardDiff.jacobian(z->c_p(z,x,param),z0)
ForwardDiff.jacobian(x->c_p(z0,x,param),x0)
ForwardDiff.jacobian(z->c_pNL(z,x,param),z0)
ForwardDiff.jacobian(x->c_pNL(z0,x,param),x0)
ForwardDiff.jacobian(z->c_(z,x,param),z0)
ForwardDiff.jacobian(x->c_(z0,x,param),x0)
ForwardDiff.jacobian(z->c_NL(z,x,param),z0)
ForwardDiff.jacobian(x->c_NL(z0,x,param),x0)

ForwardDiff.jacobian(z->c_(z,param),z0)
ForwardDiff.jacobian(z->c_NL(z,param),z0)

ForwardDiff.jacobian(z->p_NN(z,x,model),z0)
ForwardDiff.jacobian(x->p_NN(z0,x,model),x0)


ForwardDiff.gradient(z->obj(z,param),z0)

obj(z,x,op::optimParam; use_nlsolve::Bool=false)

ForwardDiff.gradient(z->obj(z,x,param),z0)
ForwardDiff.gradient(x->obj(z0,x,param),x0)

# Just do a check if eltype(u) <: AbstractFloat and only mutate if true.


# Test automatic differentiation

function f!(fvec, x)
    fvec[1] = (x[1]+3)*(x[2]^3-7)+18
    fvec[2] = sin(x[2]*exp(x[1])-1)
end
ForwardDiff.jacobian(f!,[ -0.5; 1.4])

r = nlsolve(f!, [ -0.5; 1.4], autodiff = :forward)
@test converged(r)
@test norm(r.zero - [ 0; 1]) < 1e-8





using LinearAlgebra, OrdinaryDiffEq
using DiffEqBase: get_tmp, dualcache
function foo(du, u, (A, tmp), t)
    tmp = DiffEqBase.get_tmp(tmp, u)
    mul!(tmp, A, u)
    @. du = u + tmp
    nothing
end
chunk_size = 5
prob = ODEProblem(foo, ones(5, 5), (0., 1.0), (ones(5,5), DiffEqBase.dualcache(zeros(5,5), Val{chunk_size})))
solve(prob, TRBDF2(chunk_size=chunk_size))


using ModelingToolkit, NonlinearSolve

@variables x y z
@parameters σ ρ β

# Define a nonlinear system
eqs = [0 ~ σ*(y-x),
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
ns = NonlinearSystem(eqs, [x,y,z], [σ,ρ,β])

guess = [x => 1.0,
         y => 0.0,
         z => 0.0]

ps = [
      σ => 10.0
      ρ => 26.0
      β => 8/3
      ]

prob = NonlinearProblem(ns,guess,ps)
sol = solve(prob,NewtonRaphson())

prob = NonlinearProblem(ns,guess,ps,jac=true)
sol = solve(prob,NewtonRaphson())



using ModelingToolkit, GalacticOptim, Flux, Optim

@variables x y
@parameters a b
loss = (a - x)^2 + b * (y - x^2)^2
sys = OptimizationSystem(loss,[x,y],[a,b])

u0 = [
    x=>1.0
    y=>2.0
]
p = [
    a => 6.0
    b => 7.0
]

prob = OptimizationProblem(sys,u0,p,grad=true,hess=true)
solve(prob,Newton())
sys = modelingtoolkitize(prob)
jac = eval(ModelingToolkit.generate_jacobian(sys)[2])


## Optimization

rosenbrock(x,p) =  (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
x0 = zeros(2)
p  = [1.0,100.0]

prob = OptimizationProblem(rosenbrock,x0,p)
sys = modelingtoolkitize(prob) # symbolicitize me captain!

prob = OptimizationProblem(sys,x0,p,grad=true,hess=true)
sol = solve(prob,NelderMead())
@test sol.minimum < 1e-8

sol = solve(prob,BFGS())
@test sol.minimum < 1e-8

sol = solve(prob,Newton())
@test sol.minimum < 1e-8

