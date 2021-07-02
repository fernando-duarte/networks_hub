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
mini_batch=1
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

# function betaprod_pdf!(x,a,b,f)
#     k = prod(beta.(a,b))
#     betaprod_density!(x,a,b,f)
#     f .= f/k
# end

# function betaprod_density!(x,a,b,f)
#    f[1,:] .= 1.0
#    Threads.@threads for j in 1:size(x,2)
#        for i in 1:size(x, 1)
#            f[1, j] *=  x[i,j]^(a[i] - 1.0) * (1.0 - x[i,j])^(b[i] - 1.0)
#            # f[1, j] *=   view(x, i, j:j)[]^(a[i] - 1.0) * (1.0 - view(x, i, j:j)[])^(b[i] - 1.0)
#        end
#    end
# end

function int_nlsolve!(dx,x,zp)
    z, p = zp[1:M], zp[M+1:6N+1+2N*M+M^2]
    c, α, β = selc*z, selα*z, selβ*z
    obj_x(x) = obj_nlsolve(z,p,x)
    beta_pdf(x) = betaprod_pdf(x,α,β)   # exp.(logpdf(Product(Beta.(α,β)),x))    
    dx .= obj_x(x).*beta_pdf(x)
end

function int_cont!(dx,x,zp)
  z, p = zp[1:M], zp[M+1:6N+1+2N*M+M^2]
  c, α, β = selc*z, selα*z, selβ*z
  obj_x(x) = obj_contraction(z,p,x)
  #beta_pdf(x) = betaprod_pdf(x,α,β)   # exp.(logpdf(Product(Beta.(α,β)),x))    
  beta_pdf(x) = betaprod_pdf(x,α,β)   # exp.(logpdf(Product(Beta.(α,β)),x))    
  dx .= obj_x(x).*beta_pdf(x)
end

function int_cont(x,zp)
  z, p = zp[1:M], zp[M+1:6N+1+2N*M+M^2]
  c, α, β = selc*z, selα*z, selβ*z
  obj_x(x) = obj_contraction(z,p,x)
  beta_pdf(x) = betaprod_pdf(x,α,β)   # exp.(logpdf(Product(Beta.(α,β)),x))    
  obj_x(x).*beta_pdf(x)
end

function obj_nlsolve_int(z,p)  
    mb = size(selp,1)÷N
    prob = QuadratureProblem(int_nlsolve!,zeros(N),ones(N),vcat(z,p),batch=mb)
    sol = Quadrature.solve(prob,CubaCuhre(),reltol=1e-2,abstol=1e-2)[1]
end

function obj_cont_int(z,p)
    mb = size(selp,1)÷N
    prob = QuadratureProblem(int_cont!,zeros(N),ones(N),vcat(z,p),batch=mb)
    sol = Quadrature.solve(prob,CubaCuhre(),reltol=1e-1,abstol=1e-1)[1]
end

function obj_cont_intH(z,p)
    mb = size(selp,1)÷N
    prob = QuadratureProblem(int_cont,zeros(N),ones(N),vcat(z,p))
    sol = Quadrature.solve(prob,HCubatureJL(),reltol=1e-3,abstol=1e-3)[1] # CubatureJLp, CubatureJLh, HCubatureJL
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
x0 = ones(N,mini_batch)/10