# module NetworkUtils2

# using Convex, SCS, NLsolve
# using LinearAlgebra, Distributions, Random, SpecialFunctions,  SparseArrays, MeasureTheory
# using BenchmarkTools, Test
# using DataFrames, BSON, XLSX, JLD2, Missings
# using Dates, MonthlyDates
# # using Plots, Gadfly, Cairo, Fontconfig, LightGraphs, SimpleWeightedGraphs, GraphRecipes
# using Setfield, Parameters, Suppressor, Flatten
# using MathOptInterface
# const MOI = MathOptInterface
# using Reexport

# include("IncBetaDer.jl")
# using .IncBetaDer
# @reexport using .NetworkType
# export optimParam, c_lin, c_quad, c_chance, c_p, c_pNL, c_, c_NL, c_, c_NL, obj, obj_full, loss_d



using BSON, DataFrames, LinearAlgebra, SparseArrays, Parameters, Statistics, NLsolve, Quadrature, Cuba,  DistributionsAD, SpecialFunctions, BlockArrays, Revise
import Distributions
include("NetworkType.jl"); using .NetworkType
include("IncBetaDer.jl"); #using .IncBetaDer

## instantiate
T=Float64
const N=5
mini_batch=3
data_dict = BSON.load("data.bson")
network = netEmp(data_dict[:data],5)[1]
@unpack w, p_bar, delta, γ, a = network
M = N^2+4*N+N*mini_batch

# network constants
# z = vcat(A...,b,c,α,β,p_clear...)

# constants for linear constraints
# hᵀz=q
h₁₁ᵀ = hcat(sparse(kron(Matrix{T}(I,N,N),p_bar')))
h₁ᵀ = hcat(h₁₁ᵀ,spzeros(T,N,N),sparse(I,N,N),spzeros(T,N,N),spzeros(T,N,N),spzeros(T,N,mini_batch*N))
q₁ = a

h₂₁ᵀ = repeat(spdiagm(p_bar),1,N)
h₂ᵀ = hcat(h₂₁ᵀ,sparse(I,N,N),spzeros(T,N,N),spzeros(T,N,N),spzeros(T,N,N),spzeros(T,N,mini_batch*N))
q₂ = p_bar

hᵀ = vcat(h₁ᵀ,h₂ᵀ)
q = vcat(q₁,q₂)

# constants for quadratic constraints
# zᵀHz = 0
li = LinearIndices((N,N))
liᵀ= transpose(li)
id = li[tril!(trues(N,N), 0)]
idᵀ = liᵀ[tril!(trues(N,N), 0)]
H = sparse(vcat(id,idᵀ),vcat(idᵀ,id),T(0.5),M,M)
# add a linear term to the quadratic form
# zᵀHz + h₃ᵀz = q₃
#h₃ᵀ = spzeros(T,(1,M))
#q₃ = spzeros(T,1)

const selA  = sparse(1:N^2,1:N^2,T(1),N^2,M)
const selb  = sparse(1:N,N^2+1:N^2+N,T(1),N,M)
const selc  = sparse(1:N,N^2+N+1:N^2+2N,T(1),N,M)
const selα  = sparse(1:N,N^2+2N+1:N^2+3N,T(1),N,M)
const selβ  = sparse(1:N,N^2+3N+1:N^2+4N,T(1),N,M)
const selp  = sparse(1:mini_batch*N,N^2+4N+1:M,T(1),mini_batch*N,M)


## constraints and objective
# hᵀz=q
function c_lin(z,p)
    q, hᵀ = p[4N+2:6N+1], reshape(p[6N+2:6N+1+2N*M],2N,M)
    return hᵀ*z - q  #hᵀz=q
end

# zᵀHz = 0
function c_quad(z,p)
    H  = reshape(p[6N+2+2N*M:6N+1+2N*M+M^2],M,M)
    return dot(z,H*z) #transpose(z)*H*z
end

## chance constraints ⟺ c_chance(z)=0
# Prob(c[i]x[i]>=w[i])=delta[i]
beta_ccdf(a,b,x) = zero(x) < x < one(x) ? one(x) - IncBetaDer.beta_inc_grad(a, b, x)[1] : zero(x)
function c_chance(z,p)
    c, α, β = selc*z, selα*z, selβ*z
    w, delta = p[1:N], p[2N+1:3N]
    zeroz = zero(z[1])+0.0001
    onez = one(z[1])-0.0001
    return beta_ccdf.(α,β,max.(min.(w./c,onez),zeroz) ) .- delta
end

## clearing vector constraint
# contraction map
function contraction(z,p,x)
    mb = size(x,2)
    c, p_clear, Aᵀ = selc*z, reshape(selp*z,N,mb), transpose(reshape(selA*z,N,N))
    p_bar, γ = p[N+1:2N], p[3N+1]
    p_bar_rep = repeat(p_bar,1,mb)
    c_rep = repeat(c,1,mb)
    return reshape( min.(max.(0.0,(1+γ)*(Aᵀ*p_clear + (one(x[1]) .- x).*c_rep) - γ*p_bar_rep),p_bar_rep) ,N*mb)
end

function c_p(z,p,x)
    #mb = size(x,2)
    p_clear = selp*z
    return p_clear - contraction(z,p,x)
end

function p_contraction(z, p, x, n::Integer)
    mb = size(x,2)
    p_bar, γ = p[N+1:2N], p[3N+1]
    p_bar_rep = repeat(p_bar,mb)
    znop = z[1:N^2+4N]
    n <= 0 ? p_bar_rep : contraction(vcat(znop,p_contraction(z,p,x,n-1)),p,x)
end

function p_nlsolve(z,p,x)
    mb = size(x,2)
    p_bar = p[N+1:2N]
    p_bar_rep = repeat(p_bar,mb)
    znop = z[1:N^2+4N]
    nlsolve(p_clear -> c_p(vcat(znop,p_clear),p,x),p_bar_rep).zero
end

function _c(z,p,x)
    vcat(
        c_lin(z,p),
        c_quad(z,p),
        c_chance(z,p),
        c_p(z,p,x)
        )
end

function _c(z,p)
    vcat(
        c_lin(z,p),
        c_quad(z,p),
        c_chance(z,p)
        )
end

function obj_full(obj_num,p)
    p_bar = p[N+1:2N]
    obj_num .+ sum(p_bar)
end


function obj(z,p)
    mb = size(selp,1)÷N
    c, α, β, p_clear = selc*z, selα*z, selβ*z, reshape(selp*z,N,mb)
    Ex = α./(α + β)
    return sum(c.*Ex) - sum(mean(p_clear,dims=2))
end

function obj_contraction(z,p,x)
    mb = size(x,2)
    p_clear= reshape(p_contraction(z, p, x, 3),N,mb)
    c = selc*z
    return sum(c.*x,dims=1)-sum(p_clear,dims=1) # spillovers for network
    #return sum(c.*α./(α + β))-sum(mean(p.zero,dims=2))
end

function obj_nlsolve(z,p,x)
    mb = size(x,2)
    p_clear= reshape(p_nlsolve(z,p,x),N,mb)
    c = selc*z
    return sum(c.*x,dims=1)-sum(p_clear,dims=1) # spillovers for network
    #return sum(c.*α./(α + β))-sum(mean(p.zero,dims=2))
end


function betaprod_pdf(x,a,b)
    k = prod(beta.(a,b))
    betaprod_density(x,a,b)/k
end
function betaprod_density(x,a,b)
    prod(x.^(a .- 1.0) .* (1.0 .- x).^(b .- 1.0),dims=1)
end

function betaprod_pdf!(x,a,b,f)
    k = prod(beta.(a,b))
    betaprod_density!(x,a,b,f)/k
end

function betaprod_density!(x,a,b,f)
   f[1,:] .= 1.0
   Threads.@threads for j in 1:size(x,2)
       for i in 1:size(x, 1)
           f[1, j] *=  x[i,j]^(a[i] - 1.0) * (1.0 - x[i,j])^(b[i] - 1.0)
       end
   end
end
ff=zeros(1,mini_batch)
@testset "threaded beta prob" begin
for ntest=1:100
    xtest=rand(N,mini_batch)
    betaprod_density!(xtest,selα*z0,selβ*z0,ff)
    @test isapprox(ff ,betaprod_density(xtest,selα*z0,selβ*z0))
    @test isapprox(ff ,betaprod_density(xtest,selα*z0,selβ*z0))
end
end


function int_nlsolve(dx,x,zp)
    z, p = zp[1:M], zp[M+1:6N+1+2N*M+M^2]
    c, α, β = selc*z, selα*z, selβ*z
    obj_x(x) = obj_nlsolve(z,p,x)
    beta_pdf(x) = betaprod_pdf(x,α,β)   # exp.(logpdf(Product(Beta.(α,β)),x))
    dx .= obj_x(x).*beta_pdf(x)
end

function int_cont(dx,x,zp)
  z, p = zp[1:M], zp[M+1:6N+1+2N*M+M^2]
  c, α, β = selc*z, selα*z, selβ*z
  obj_x(x) = obj_contraction(z,p,x)
  #beta_pdf(x) = betaprod_pdf(x,α,β)   # exp.(logpdf(Product(Beta.(α,β)),x))
  beta_pdf(x) = betaprod_pdf(x,α,β)   # exp.(logpdf(Product(Beta.(α,β)),x))
  dx .= obj_x(x).*beta_pdf(x)
end

function obj_nlsolve_int(z,p)
    mb = size(selp,1)÷N
    prob = QuadratureProblem(int_nlsolve,zeros(N),ones(N),vcat(z,p),batch=mb)
    sol = Quadrature.solve(prob,CubaCuhre(),reltol=1e-2,abstol=1e-2)[1]
end


function obj_cont_int(z,p)
    mb = size(selp,1)÷N
    prob = QuadratureProblem(int_cont,zeros(N),ones(N),vcat(z,p),batch=mb)
    sol = Quadrature.solve(prob,CubaCuhre(),reltol=1e-1,abstol=1e-1)[1]
end

# upper and lower bounds
# low<=z<=upp
lowα = 0.5 ;uppα = 10
lowβ = 1   ;uppβ = 30
low_zero = eps(T) # eps(T) T(0.005)
low = vcat(spzeros(T,N^2),fill(low_zero,N),fill(low_zero,N),fill(lowα,N),fill(lowβ,N),fill(low_zero,N*mini_batch))
low = Array(low)
upp = vcat(ones(T,N^2),p_bar,a,Array{T}(fill(uppα,N)),Array{T}(fill(uppβ,N)),repeat(p_bar,mini_batch))

p_sparse = vcat(w, p_bar, delta, γ, a ,q ,hᵀ...,H...)
p_dense =  vcat(w, p_bar, delta, γ, a ,q,Array(hᵀ)...,Array(H)...)

z0 = low + (upp-low)/2;
A0 = (LowerTriangular(ones(N,N))-Matrix(I,N,N))/2
z0[1:N^2] = [A0...]
p0 = p_dense

c_lin(z0,p0)
c_quad(z0,p0)
c_chance(z0,p0)
_c(z0,p0)

x0 = ones(N,mini_batch)/10
c_p(z0,p0,x0)
_c(z0,p0,x0)

contraction(z0,p0,x0)
p1=p_contraction(z0, p0, x0, 10)
p2=p_nlsolve(z0,p0,x0)
p1 ≈ p2

obj(z0,p0)
obj_full(obj(z0,p0),p0)
obj_full(obj_contraction(z0,p0,x0),p0)
obj_full(obj_nlsolve(z0,p0,x0),p0)


obj_contraction(z0,p0,x0)
obj_nlsolve(z0,p0,x0)
α,β= selα*z0,selβ*z0
betaprod_pdf(x0,α,β)
betaprod_density(x0,α,β)
dx0=zeros(1,mini_batch)
int_nlsolve(dx0,x0,vcat(z0,p0))
int_cont(dx0,x0,vcat(z0,p0))
obj_nlsolve_int(z0,p0)
obj_cont_int(z0,p0)


# integrate by drawing x
dist_uni = Distributions.Product(Distributions.Uniform.(zeros(size(α)),ones(size(β))))
dist_beta = Distributions.Product(Distributions.Beta.(α,β))
num_draw = 100000
num_draw = (num_draw ÷ mini_batch)*mini_batch # a multiple of mini_batch

uni_draw = rand(dist_uni,num_draw)
beta_draw = rand(dist_beta,num_draw)
@show isapprox(mean(uni_draw,dims=2),fill(0.5,N), atol=1e-3);
@show isapprox(mean(beta_draw,dims=2),α./(α + β), atol=1e-3);

uni_block = BlockArray(uni_draw,[N],fill(mini_batch,num_draw ÷ mini_batch))
beta_block = BlockArray(beta_draw,[N],fill(mini_batch,num_draw ÷ mini_batch))

# using contraction mapping
int_cont0(x) = int_cont(dx0,x,vcat(z0,p0))
obj_cont0(x) = obj_contraction(z0,p0,x)
approx_uni  = mean([mean(int_cont0(view(uni_block, Block(i)))) for i=1:blocklength(uni_block)])
approx_beta = mean([mean(obj_cont0(view(beta_block, Block(i)))) for i=1:blocklength(beta_block)])
approx_int  = obj_cont_int(z0,p0)
@show approx_uni, approx_beta, approx_int ;


# using nlsolve
int_nlsolve0(x) = int_nlsolve(dx0,x,vcat(z0,p0))
obj_nlsolve0(x) = obj_nlsolve(z0,p0,x)
approx_uni = mean([mean(int_nlsolve0(view(uni_block, Block(i)))) for i=1:blocklength(uni_block)])
approx_beta =mean([mean(obj_nlsolve0(view(beta_block, Block(i)))) for i=1:blocklength(beta_block)])
approx_int = obj_nlsolve_int(z0,p0)
@show approx_uni, approx_beta, approx_int ;

## JuMP
using JuMP, Ipopt
m = Model(with_optimizer(Ipopt.Optimizer,
                            max_iter=100000,
                            print_level=5,
                            mu_strategy="adaptive",
                            derivative_test="second-order",
                            hessian_approximation="exact"
                        )
         )

# @variable(m, low[i] <= z[i=1:M] <= upp[i], start = z0[i])
# @objective(m, Max,sum(z[i]^2 for i=1:M))
# JuMP.optimize!(m)

## constraints and objective
# hᵀz=q
qJ, hᵀJ = p0[4N+2:6N+1], reshape(p0[6N+2:6N+1+2N*M],2N,M)
@constraint(m, hᵀJ*z .== qJ )

# zᵀHz = 0
HJ  = reshape(p0[6N+2+2N*M:6N+1+2N*M+M^2],M,M)
@constraint(m,dot(z,HJ*z) == 0 )

for i=1:N
    @eval ($(Symbol("c_chanceJ$i")))(z...) = c_chance(collect(z),p0)[$i]
    @eval JuMP.register(m, :($(Symbol("c_chanceJ$i"))) , M,  ($(Symbol("c_chanceJ$i"))), autodiff=true)
    @eval @NLconstraint(m,($(Symbol("c_chanceJ$i")))(z...) == 0)
end

objJ(z...)=obj_cont_int(collect(z),p0)
JuMP.register(m, :objJ, M, objJ, autodiff=true)
@NLobjective(m, Min , objJ(z...)  )

unset_silent(m)
JuMP.optimize!(m)

st0 = termination_status(m)
obj0 = objective_value(m)
@show st0, obj0

Asol0 = reshape(selA*JuMP.value.(z),N,N)
bsol0 = selb*JuMP.value.(z)
csol0 = selc*JuMP.value.(z)
αsol0 = selα*JuMP.value.(z)
βsol0 = selβ*JuMP.value.(z)
psol0 = reshape(selp*JuMP.value.(z),N,mini_batch)

zsol0 = [Asol0...;bsol0...;csol0...;psol0...]
Zsol0 = zsol0*zsol0'
rank(Zsol0)
eigvals(Zsol0)


cdf_pkg = [1-cdf(Beta(αsol0[i],βsol0[i]),w[i]/csol0[i]) for i=1:N]
cdf_opt = [dist_cdf([αsol0[i],βsol0[i],w[i]/csol0[i] ]...) for i=1:N]
tol = 1e-5
@testset "check solution" begin
    @test norm( sum(Asol0,dims=2).* data.p_bar .- (data.p_bar .- bsol0)) < tol
    @test norm( Asol0' * data.p_bar .- (data.assets .- csol0)) < tol
    @test norm(diag(Asol0)) < tol
    @test norm([Asol0[i,j]*Asol0[j,i] for i=1:N , j=1:N]) < tol
    @test all(0 .<=Asol0.<=1)
    @test all(0 .<=bsol0.<=data.p_bar)
    @test all(0 .<=csol0.<=data.assets)
    @test all(0 .<=psol0.<=data.p_bar)

    @test norm(min.(p_bar,max.((1+g0).*(Asol0'*psol0[:,1] .+ csol0 .- x0[:,1].*csol0 ) .- g0.*p_bar,0))-psol0[:,1])<tol


    @test norm(cdf_pkg-cdf_opt)<0.01
    @test norm(delta-cdf_opt)<0.01

end


## NLP Models
using NLPModels, ADNLPModels, NLPModelsIpopt, ReverseDiff, NLPModelsModifiers
nc = length(_c(z0,p0) )
nlp = ADNLPModel(z->obj_cont_int(z,p0), z0, low, upp, z->_c(z,p0), zeros(nc), zeros(nc); lin=1:2N) #  adbackend = ADNLPModels.ReverseDiffAD()  adbackend = ADNLPModels.ZygoteAD() ; lin=1:2N)
obj(nlp,nlp.meta.x0)
gx = grad(nlp, nlp.meta.x0)
nlpLBFGS=NLPModelsModifiers.LBFGSModel(nlp)

function hess(nlp::AbstractNLPModel, x::AbstractVector; obj_weight::Real = one(eltype(x)))
  @lencheck nlp.meta.nvar x
  rows, cols = hess_structure(nlp)
  vals = hess_coord(nlp, x, obj_weight = obj_weight)
  Symmetric(sparse(FiniteDiff.finite_difference_hessian(x->obj(nlp,x),x)), :L)
end

hess(nlpLBFGS,z0)


Hx = hess(nlp, nlp.meta.x0)

obj(nlp, nlp.meta.x0) # evaluate the objective value at x0
grad!(nlp, nlp.meta.x0, g) # evaluate the objective gradient at x0
cons!(nlp, nlp.meta.x0, c) # evaluate the vector of constraints at x0


jac_structure!(nlp, rows, cols): fill rows and cols with the spartity structure of the Jacobian, if the problem is constrained;
jac_coord!(nlp, nlp.meta.x0, vals): fill vals with the Jacobian values corresponding to the sparsity structure returned by jac_structure!();
hess_structure!(nlp, rows, cols): fill rows and cols with the spartity structure of the lower triangle of the Hessian of the Lagrangian;
hess_coord!(nlp, x, y, vals; obj_weight=1.0): fill vals with the values of the Hessian of the Lagrangian corresponding to the sparsity structure returned by hess_structure!(), where obj_weight is the weight assigned to the objective, and y is the vector of multipliers.

stats = ipopt(nlpLBFGS)
print(stats)
problem = createProblem(n, nlp.meta.lvar, nlp.meta.uvar,
m, nlp.meta.lcon, nlp.meta.ucon,
nlp.meta.nnzj, nlp.meta.nnzh,
eval_f, eval_g, eval_grad_f,
eval_jac_g);#, eval_h)



using DCISolver, ADNLPModels, Logging


grad(nlp,z0) == ForwardDiff.gradient(z->obj_cont_int(z,p0),z0)
dci(nlp, nlp.meta.x0)

stats = with_logger(NullLogger()) do
  dci(nlp, nlp.meta.x0)
end

println(stats)




Pkg.add(url="https://github.com/JuliaSmoothOptimizers/NCL.jl", rev="master")
using NLPModelsIpopt, NCL
nlp = ADNLPModel(z->obj_cont_int(z,p0), z0, low, upp, z->_c(z,p0), zeros(nc), ones(nc)) # adbackend = ADNLPModels.ZygoteAD() ; lin=1:2N)

ncl_nlin_res = NCLModel(nlp)
NCLSolve(ncl_nlin_res)








# end #module



using PolyChaos
N=5
α=ones(N);β=2.0*ones(N);
mop = MultiOrthoPoly([ Beta01OrthoPoly(5, α[i],β[i]) for i in 1:N])







import Pkg; Pkg.add("SparseGrids");Pkg.add("FastGaussQuadrature");Pkg.add("PolyChaos")
using SparseGrids, FastGaussQuadrature, PolyChaos
#If f is a function that returns nodes, weights = f(n), for any (integer) order n, then the function sparsegrid computes the sparse extension to D dimensions of order O
nodes, weigths = sparsegrid(D, O, f)

order = [3, 5, 6, 4]
op1 = GaussOrthoPoly(order[1], addQuadrature=true);
op2 = Uniform01OrthoPoly(order[2], addQuadrature=true);
op3 = Beta01OrthoPoly(order[3], 2, 1.2, addQuadrature=true);
ops = [op1, op2, op3];
integrate(c->integrate(b->integrate((a)->a*b*c, op1),op2),op3)


order, Nrec, N, tol = 2, 3, 3, 1e-4
mop = MultiOrthoPoly([ Beta01OrthoPoly(order, α[i],β[i];  Nrec=Nrec) for i in 1:N], order)
α =rand(3) .+ 1.0; β=rand(3).+2.0
op1=Beta01OrthoPoly(order, α[1],β[1];  Nrec=Nrec)
op2=Beta01OrthoPoly(order, α[2],β[2];  Nrec=Nrec)
op2=Beta01OrthoPoly(order, α[3],β[3];  Nrec=Nrec)
function f()
    integrate(c->integrate(b->integrate((a)->a, op1),op2),op3)
end
α[1]/(α[1]+β[1])
@btime $f()

show(mop)
calculateAffinePCE(mop)
Quad(10,mop.measure)
a,b = coeffs(mop)
mean(coeffs(mop),mop)
calculateAffinePCE(1,mop)

# all-purpose constructor (last resort!)
function Quad(N::Int, w::Function, dom::Tuple{<:Real,<:Real}; quadrature::Function=clenshaw_curtis)
    N <= 0 && throw(DomainError(N, "number of quadrature points has to be positive"))
    nodes, weights = quadgp(w, dom[1], dom[2], N; quadrature=quadrature)
    Quad("quadgp", N, nodes, weights)
end

function Quad(N::Int, measure::AbstractMeasure; quadrature::Function=clenshaw_curtis)
    typeof(measure) != Measure && @warn "For measures of type $(typeof(measure)) the quadrature rule should be based on the recurrence coefficients."
    Quad(N, measure.w, (measure.dom[1], measure.dom[2]); quadrature=quadrature)
end







Tensor(2,mop)

PolyChaos.integrate(x-> x[1]+x[2], mop)

nodes, weights = sparsegrid(3,20,n->nw(n,1.0,2.0))
# Evaluate integrand
F = map(x -> x, nodes)
Fn = map(sum, F)

# Evaluate quadrature sum
Q = dot(Fn, weights)



using PolyChaos, Test, LinearAlgebra
import QuadGK: quadgk

function integration_test(op::AbstractOrthoPoly, dim::Int,name::String;tol::Float64=1e-6,
                          dom::Tuple{<:Real,<:Real}=op.measure.dom)
    N, w, ind = op.deg, op.measure.w, zeros(Int64,dim)
    @testset "$name" begin
        for ind_ in Iterators.product([collect(1:N) for i in 1:dim]...)
            ind[:] .= ind_[:]
            s1 = computeSP(ind, op)
            f(t) = prod(evaluate(ind[i], t, op) for i in 1:dim)*w(t)
            s2 = quadgk(f, dom[1], dom[2])[1]
            if abs(s1)<=1e-3
                @test isapprox(s1, s2; atol=tol)
            else
                @test isapprox(abs(s1/s2), 1; atol=tol)
            end
        end
    end
end

N, Nrec, dim, tol = 3, 1000, 3, 1e-4
α, β = 1.32, 3.23
op = Beta01OrthoPoly(N, α, β; Nrec=Nrec)
set_beta = integration_test(op, dim, "beta01"; tol=tol)
N, w, ind = op.deg, op.measure.w, zeros(Int64,dim)
ind_ = (1,2,3)
dom = (0,1)

for ind_ in Iterators.product([collect(1:N) for i in 1:dim]...)
    ind[:] .= ind_[:]
    s1 = computeSP(ind, op)
    f(t) = prod(evaluate(ind[i], t, op) for i in 1:dim)*w(t)
    s2 = quadgk(f, dom[1], dom[2])[1]
    if abs(s1)<=1e-3
        @test isapprox(s1, s2; atol=tol)
    else
        @test isapprox(abs(s1/s2), 1; atol=tol)
    end
end

for ind_ in Iterators.product([collect(1:N) for i in 1:dim]...)
   display(ind_[:] )
end



using PolyChaos
dd = [3]
op1 = JacobiOrthoPoly(dd[1],1.1,2.3)
op2 = JacobiOrthoPoly(dd[1],1.1,2.3)
show(op1)
ops = [op1,op2]
mop = MultiOrthoPoly(ops, dd[1])
mop.ind

using Revise
using SmolyakApprox, FastGaussQuadrature, LinearAlgebra
d  = 1  # Set the number of dimensions
mu = 3  # Set the level of approximation
g1, m1 = smolyak_grid(chebyshev_gauss_lobatto,d,mu)  # Construct the Smolyak grid and the multi index
g2, m2 =  smolyak_grid(clenshaw_curtis_equidistant,d,mu)
gj(n::Integer) = gaussjacobi(n,1.1,2.3)
g3, m3 = smolyak_grid(gj,d,mu)

function test(grid)
    #y_value = (grid[:,1].+1).^0.1.*exp.(grid[:,2]).*log.(grid[:,3].+2).^0.2.*(grid[:,4].+2).^0.8.*(grid[:,5].+7).^0.1
    #y_value = (grid[:,1].+1).^0.1.*exp.(grid[:,2]).*log.(grid[:,1].+2).^0.2.*(grid[:,2].+2).^0.8.*(grid[:,1].+7).^0.1
    y_value = (grid[:,1].+1).^0.1
    return y_value
end
y1 = test(g1)
y2 = test(g2)
y3 = test(g3)

w1 = smolyak_weights(y1,g1,m1)
w2 = smolyak_weights(y2,g2,m2)
w3 = smolyak_weights(y3,g3,m3)


im1 = smolyak_inverse_interpolation_matrix(g1,m1)  # Compute the interpolation matrix
im2 = smolyak_inverse_interpolation_matrix(g2,m2)  # Compute the interpolation matrix
im3 = smolyak_inverse_interpolation_matrix(g3,m3)  # Compute the interpolation matrix

w1 = smolyak_weights(y1,im1)             # Compute the Smolyak weights
w2 = smolyak_weights(y2,im2)             # Compute the Smolyak weights
w3 = smolyak_weights(y3,im3)             # Compute the Smolyak weights

#point = [0.75, 0.45, 0.82, -0.15, -0.95]
point = [0.75]
y_hatt1 = smolyak_evaluate(w1,point,m1)  # Evaluate the approximated function
y_hatt2 = smolyak_evaluate(w2,point,m2)  # Evaluate the approximated function
y_hatt3 = smolyak_evaluate(w3,point,m3)  # Evaluate the approximated function

  multi_index        = SmolyakApprox.generate_multi_index(d,mu)
  unique_multi_index = sort(unique(multi_index))
  unique_node_number = SmolyakApprox.m_i(unique_multi_index)

T= Float64
  base_nodes   = Array{Array{T,1},1}(undef,length(unique_node_number))
  base_weights = Array{Array{T,1},1}(undef,length(unique_node_number))
  base_nodes_gj   = Array{Array{T,1},1}(undef,length(unique_node_number))
  base_weights_gj = Array{Array{T,1},1}(undef,length(unique_node_number))
  for i = 1:length(unique_node_number)
    base_nodes_gj[i], base_weights_gj[i] = gj(unique_node_number[i])
    base_nodes[i], base_weights[i] = chebyshev_gauss_lobatto(unique_node_number[i])
  end

  unique_base_nodes_gj = Array{Array{T,1},1}(undef,length(unique_node_number))
  unique_base_nodes_gj[1] = base_nodes_gj[1]
  unique_base_nodes = Array{Array{T,1},1}(undef,length(unique_node_number))
  unique_base_nodes[1] = base_nodes[1]
  for i = 2:length(unique_base_nodes)
    unique_base_nodes[i] = setdiff(base_nodes[i],base_nodes[i-1])
    unique_base_nodes_gj[i] = setdiff(base_nodes_gj[i],base_nodes_gj[i-1])
  end


  nodes = Array{T,2}(undef,determine_grid_size(multi_index))
  #nodes_gj = Array{T,2}(undef,determine_grid_size(multi_index))
  nodes_gj = Array{T,2}(undef,(sum(map(x->length(x),unique_base_nodes_gj)), d))
  l = 1
  l_gj = 1
  @inbounds for j = 1:size(multi_index,1)
    new_nodes = unique_base_nodes[multi_index[j,1]]  # Here new_nodes is a 1d array
    new_nodes_gj = unique_base_nodes_gj[multi_index[j,1]]  # Here new_nodes is a 1d array
    for i = 2:d
      new_nodes = combine_nodes(new_nodes,unique_base_nodes[multi_index[j,i]])  # Here new_nodes becomes a 2d array
      new_nodes_gj = combine_nodes(new_nodes_gj,unique_base_nodes_gj[multi_index[j,i]])  # Here new_nodes becomes a 2d array
    end
    m = size(new_nodes,1)
    m_gj = size(new_nodes_gj,1)
    nodes[l:l+m-1,:] = new_nodes
    nodes_gj[l_gj:l_gj+m_gj-1,:] = new_nodes_gj
    l += m
    l_gj += m_gj
  end






function determine_grid_size(mi)

  temp = similar(mi)

  for i = 1:size(mi,1)
    for j = 1:size(mi,2)
      if mi[i,j] == 1
        temp[i,j] = 1
      elseif mi[i,j] == 2
        temp[i,j] = 2^(mi[i,j]-1)
      else
        temp[i,j] = 2^(mi[i,j]-1)+1 - (2^(mi[i,j]-2)+1)
      end
    end
  end

  s = 0
  for i = 1:size(mi,1)
    t = 1
    for j = 1:size(mi,2)
      t *= temp[i,j]
    end
    s += t
  end

  return (s, size(mi,2))

end
#= A second way of computing the weights and evaluating the approximated function that
  computes the interpolation matrix just once. =#

  interp_mat = smolyak_inverse_interpolation_matrix(grid,multi_ind)  # Compute the interpolation matrix
  w = smolyak_weights(y,interp_mat)             # Compute the Smolyak weights
  y_hatt = smolyak_evaluate(w,point,multi_ind)  # Evaluate the approximated function

  # Piecewise linear
 y_pl = test(g2)
  w_pl = smolyak_pl_weights(y_pl,g2,m2)
  y_pl_hat = smolyak_pl_evaluate(w_pl,point,g2,m2)

  # Evaluate the exact function at point

  y_actual = test(point')
y_hatt


function nw(n::S,domain = [1.0,-1.0])  where {S<:Integer}
    a=1.0
    b=2.0
    opq = Beta01OrthoPoly(n, a, b, addQuadrature = true)
    opq = Beta01OrthoPoly(n, a, b, addQuadrature = true)
    jacobi
    return (opq.quad.nodes,opq.quad.weights)
end

clenshaw_curtis_equidistant(5)
nw(5)

y = grid[:,1]
w = smolyak_weights(y,interp_mat)

interp_mat = smolyak_inverse_interpolation_matrix(grid,multi_ind)  # Compute the interpolation matrix


function chebyshev_gauss_lobatto(n::S,domain = [1.0,-1.0]) where {S<:Integer}

  # These nodes a just the Chebyshev extrema by a different name.

  # Construct the nodes on the [-1.0,1.0] interval

  if n <= 0
    error("The number of nodes must be positive.")
  end

  if n == 1
    nodes   = [0.0]
    weights = [1.0*pi]
  else
    nodes    = zeros(n)
    nodes[1] = -1.0
    nodes[n] = 1.0

    weights    = zeros(n)
    weights[1] = pi/(2*(n-1))
    weights[n] = pi/(2*(n-1))

    for i = 2:div(n,2)
      x = cos(pi*(i-1)/(n-1))
      nodes[i]       = -x
      nodes[end-i+1] = x
      weights[i]       = pi/(n-1)
      weights[end-i+1] = pi/(n-1)
    end

    if isodd(n)
      nodes[round(Int,(n+1)/2)]   = 0.0
      weights[round(Int,(n+1)/2)] = pi/(n-1)
    end
  end

  # Scale the nodes to the desired domain

  nodes = scale_nodes(nodes,domain)

  return nodes, weights

end







using Surrogates, PolyChaos
using Plots
n = 20
lower_bound = 1.0
upper_bound = 6.0
x = sample(n,lower_bound,upper_bound,LowDiscrepancySample(2))
f = x -> log(x)*x + sin(x)
y = f.(x)
scatter(x, y, label="Sampled points", xlims=(lower_bound, upper_bound), legend=:top)
plot!(f, label="True function", xlims=(lower_bound, upper_bound), legend=:top)

poly1 = PolynomialChaosSurrogate(x,y,lower_bound,upper_bound)
poly2 = PolynomialChaosSurrogate(x,y,lower_bound,upper_bound, op = GaussOrthoPoly(5))
plot(x, y, seriestype=:scatter, label="Sampled points", xlims=(lower_bound, upper_bound), legend=:top)
plot!(f, label="True function",  xlims=(lower_bound, upper_bound), legend=:top)
plot!(poly1, label="First polynomial",  xlims=(lower_bound, upper_bound), legend=:top)
plot!(poly2, label="Second polynomial",  xlims=(lower_bound, upper_bound), legend=:top)




## test quadrature
using Quadrature, ForwardDiff, FiniteDiff, Zygote, Cuba
f(x,p) = sum(sin.(x .* p))
lb = ones(2)
ub = 3ones(2)
p = [1.5,2.0]

function testf(p)
    prob = QuadratureProblem(f,lb,ub,p)
    sin(solve(prob,CubaCuhre(),reltol=1e-6,abstol=1e-6)[1])
end
dp1 = Zygote.gradient(testf,p)
dp2 = FiniteDiff.finite_difference_gradient(testf,p)
dp3 = ForwardDiff.gradient(testf,p)
dp1[1] ≈ dp2 ≈ dp3

hp1 = Zygote.hessian(testf,p)
hp2 = FiniteDiff.finite_difference_hessian(testf,p)
hp3 = ForwardDiff.hessian(testf,p)

Zygote.jacobian(p->Zygote.gradient(testf,p),p)
 jacobian(x -> gradient(testf, x)[1], p)[1]
ForwardDiff.jacobian(x->ForwardDiff.gradient(testf,x),p)

hess1 = ForwardDiff.jacobian(x->FiniteDiff.finite_difference_gradient(testf,x),p)
hess2 = FiniteDiff.finite_difference_jacobian(x->ForwardDiff.gradient(testf,x),p)
hess3 = FiniteDiff.finite_difference_hessian(testf,p)