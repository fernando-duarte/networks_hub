module NetDefs

using DataFrames, LinearAlgebra, SparseArrays, Parameters, Statistics, NLsolve, Quadrature, Cuba, SpecialFunctions, Distributions
include("NetworkType.jl"); using .NetworkType
#include("IncBetaDer.jl"); #using .IncBetaDer
include("ExtendedPowerDist.jl");
using SparseGrids, FastGaussQuadrature

using Flux, BSON
using BSON: @load, @save
using Flux.Data: DataLoader;using IterTools: ncycle;using Flux: throttle

export T, N, M, P, p0, z0, x0, c_lin,c_quad,c_chance, c_p, low, upp, cL, uL, contraction, p_contraction, obj_contraction, obj_contraction_uniform, mini_batch, obj_contraction_quaduniform, weights0, sump_contraction, extpow_prod_pdf, selA, selb, selc, selα, selβ

## instantiate
T=Float64
const N=5
# quadrature for uniform distribution
const nodes0, weights0 = sparsegrid(N, 2, gausslegendre)
x0 = (hcat(nodes0...) .+ 1) ./ 2 # rescale quadrature nodes from [-1,1] to [0,1]
x0 = ones(size(x0))/10
mini_batch=size(x0,2)
data_dict = BSON.load("data.bson")
network = netEmp(data_dict[:data],N)[1]
@unpack w, p_bar, delta, γ, a = network
M = N^2+4*N

#=
BSON.load(joinpath(@__DIR__,"clearing_p_NN_gpu_sum_N$N.bson"), @__MODULE__) #loads NN called `model`
#model(vcat(z0[1:N^2+2*N],x0[:,1])) # #model(vcat(A...,b,c,x)) gives sum of elements of clearing vector p

function obj_NN(z,p,x)
    mb = size(x,2)
    zNN = z[1:N^2+2*N]
    sump_clear= transpose([model(vcat(zNN,x[:,i]))[1] for i=1:mb])
    c = selc*z
    return sum(c.*x,dims=1) - sump_clear # spillovers for network
    #return sum(c.*α./(α + β))-sum(mean(p.zero,dims=2))
end
=#

# network constants
# z = vcat(A...,b,c,α,β)

# constants for linear constraints
# hᵀz=q
const h₁₁ᵀ = hcat(sparse(kron(Matrix{T}(I,N,N),p_bar')))
const h₁ᵀ = hcat(h₁₁ᵀ,spzeros(T,N,N),sparse(I,N,N),spzeros(T,N,N),spzeros(T,N,N))
const q₁ = a

const h₂₁ᵀ = repeat(spdiagm(p_bar),1,N)
const h₂ᵀ = hcat(h₂₁ᵀ,sparse(I,N,N),spzeros(T,N,N),spzeros(T,N,N),spzeros(T,N,N))
const q₂ = p_bar

const hᵀ = vcat(h₁ᵀ,h₂ᵀ)
const q = vcat(q₁,q₂)

# constants for quadratic constraints
# zᵀHz = 0
li = LinearIndices((N,N))
liᵀ= transpose(li)
id = li[tril!(trues(N,N), 0)]
idᵀ = liᵀ[tril!(trues(N,N), 0)]
const H = sparse(vcat(id,idᵀ),vcat(idᵀ,id),T(0.5),M,M)
# add a linear term to the quadratic form
# zᵀHz + h₃ᵀz = q₃
#h₃ᵀ = spzeros(T,(1,M))
#q₃ = spzeros(T,1)

# TO DO: replace by indexing?
const selA  = sparse(1:N^2,1:N^2,T(1),N^2,M)
const selb  = sparse(1:N,N^2+1:N^2+N,T(1),N,M)
const selc  = sparse(1:N,N^2+N+1:N^2+2N,T(1),N,M)
const selα  = sparse(1:N,N^2+2N+1:N^2+3N,T(1),N,M)
const selβ  = sparse(1:N,N^2+3N+1:N^2+4N,T(1),N,M)

# dense 
# const selAd  = Array(selA)
# const selbd  = Array(selb)
# const selcd  = Array(selc)
# const selαd  = Array(selα)
# const selβd  = Array(selβ)

## constraints and objective
# hᵀz=q
function c_lin(z,p)
    q, hᵀ = p[4N+2:6N+1], reshape(p[6N+2:6N+1+2N*M],2N,M)
    return hᵀ*z - q  #hᵀz=q
end

# zᵀHz = 0
function c_quad(z,p)
    H  = reshape(p[6N+2+2N*M:6N+1+2N*M+M^2],M,M)
    return dot(z,H*z) # try dot(z, H, z) #transpose(z)*H*z
end

## chance constraints ⟺ c_chance(z)=0
# Prob(c[i]x[i]>=w[i])=delta[i]
dist_ccdf(a,b,x) = 1 - ExtendedPowerDist.cdf(ExtendedPowerDist.ExtPow(a,b),x)

function c_chance(z,p)
    c, α, β = selc*z, selα*z, selβ*z
    w, delta = p[1:N], p[2N+1:3N]
    return dist_ccdf.(α,β,max.(min.(w./c,1),0)) .- delta
end

## clearing vector constraint
# contraction map
# function contraction(z,p,x)
#     mb = size(x,2)
#     c, p_clear, Aᵀ = selc*z, reshape(selp*z,N,mb), transpose(reshape(selA*z,N,N))
#     p_bar, γ = p[N+1:2N], p[3N+1]
#     p_bar_rep = repeat(p_bar,1,mb)
#     c_rep = repeat(c,1,mb)
#     return reshape( min.( (1+γ)*(Aᵀ*p_clear + (1 .- x).*c_rep) - γ*p_bar_rep, p_bar_rep) ,N*mb)
# end

function p_contraction(z,p,x)
    mb = size(x,2)
    c, Aᵀ = selc*z, transpose(reshape(selA*z,N,N))
    p_bar, γ = p[N+1:2N], p[3N+1]
    p_bar_rep = repeat(p_bar,1,mb)
    c_rep = repeat(c,1,mb)
    p1 = min.( (1+γ)*(Aᵀ*p_bar_rep + (1 .- x).*c_rep) - γ*p_bar_rep, p_bar_rep)
    #p2 = min.( (1+γ)*(Aᵀ*p1        + (1 .- x).*c_rep) - γ*p_bar_rep, p_bar_rep)
    #p3 = min.( (1+γ)*(Aᵀ*p2        + (1 .- x).*c_rep) - γ*p_bar_rep, p_bar_rep)
    #p_n = min.( (1+γ)*(Aᵀ*p_{n-1} + (1 .- x).*c_rep) - γ*p_bar_rep, p_bar_rep)
    return p1
end

function sump_contraction(z,p,x)
    mb = size(x,2)
    c, Aᵀ = selc*z, transpose(reshape(selA*z,N,N))
    p_bar, γ = p[N+1:2N], p[3N+1]
    p_bar_rep = repeat(p_bar,1,mb)
    c_rep = repeat(c,1,mb)
    p1 = min.( (1+γ)*(Aᵀ*p_bar_rep + (1 .- x).*c_rep) - γ*p_bar_rep, p_bar_rep)
    p2 = min.( (1+γ)*(Aᵀ*p1        + (1 .- x).*c_rep) - γ*p_bar_rep, p_bar_rep)
    p3 = min.( (1+γ)*(Aᵀ*p2        + (1 .- x).*c_rep) - γ*p_bar_rep, p_bar_rep)
    #p_n = min.( (1+γ)*(Aᵀ*p_{n-1} + (1 .- x).*c_rep) - γ*p_bar_rep, p_bar_rep)
    return sum(p3,dims=1)
end


# function c_p(z,p,x)
#     #mb = size(x,2)
#     p_clear = selp*z
#     return p_clear - contraction(z,p,x)
# end


function p_contraction(z, p, x, n::Integer)
    mb = size(x,2)
    p_bar, γ = p[N+1:2N], p[3N+1]
    p_bar_rep = repeat(p_bar,mb)
    znop = z[1:N^2+4N]
    n <= 0 ? p_bar_rep : contraction(vcat(znop,p_contraction(z,p,x,n-1)),p,x)
end

# function p_nlsolve(z,p,x)
#     mb = size(x,2)
#     p_bar = p[N+1:2N]
#     p_bar_rep = repeat(p_bar,mb)
#     znop = z[1:N^2+4N]
#     nlsolve(p_clear -> c_p(vcat(znop,p_clear),p,x),p_bar_rep).zero
# end

function _c(z,p,x)
    vcat(
        c_lin(z,p),
        c_quad(z,p),
        c_p(z,p,x),
        c_chance(z,p)
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

# function obj(z,p)
#     mb = size(selp,1)÷N
#     c, α, β, p_clear = selc*z, selα*z, selβ*z, reshape(selp*z,N,mb)
#     Ex = α./(α + β)
#     return sum(c.*Ex) - sum(mean(p_clear,dims=2))
# end

function ev(α,β)
   return 1 .- (1/2)*sqrt.(pi ./ β).*exp.( (α .+ 1).^2 ./ (4β)).*erfc.( (α .+ 1) ./ (2sqrt.(β)) ) 
end

function obj_contraction(z,p,x)
    # x drawn from extended power dist / beta dist
    c, α, β = selc*z, selα*z, selβ*z
    return sum(c.*ev(α,β)) - mean(sump_contraction(z,p,x)) # spillovers for network
end

# function betaprod_pdf(x,a,b)
#     k = prod(beta.(a,b))
#     betaprod_density(x,a,b)/k
# end
# function betaprod_density(x,a,b)
#     prod(x.^(a .- 1.0) .* (1.0 .- x).^(b .- 1.0),dims=1)
# end


function extpow_prod_pdf(x,a,b)
    prod(ExtendedPowerDist.density.(ExtendedPowerDist.ExtPow.(a,b),x),dims=1)
end
function extpow_prod_pdf(x,z)
    α, β = selα*z, selβ*z
    extpow_prod_pdf(x,α,β)
end

function obj_contraction_uniform(z,p,x)
    # x drawn from uniform distribution
    c, α, β = selc*z, selα*z, selβ*z
    return sum(c.*ev(α,β)) - mean(sump_contraction(z,p,x).*extpow_prod_pdf(x,α,β)) # spillovers for network
end

function obj_contraction_quaduniform(z,p,x,weights)
    # quadrature for uniform distribution (x are nodes)
    c, α, β = selc*z, selα*z, selβ*z
    return sum(c.*ev(α,β)) - dot(sump_contraction(z,p,x).*extpow_prod_pdf(x,α,β), weights)/2^N # spillovers for network
end

function obj_contraction_quad(z,p,x,weigths)
     # quadrature for extended power distribution (x are nodes)
    c, α, β = selc*z, selα*z, selβ*z
    return sum(c.*ev(α,β)) - dot(sump_contraction(z,p,x),weigths) # spillovers for network
end

function obj_nlsolve(z,p,x)
    mb = size(x,2)
    p_clear= reshape(p_nlsolve(z,p,x),N,mb)
    c = selc*z
    return sum(c.*x,dims=1)-sum(p_clear,dims=1) # spillovers for network
    #return sum(c.*α./(α + β))-sum(mean(p.zero,dims=2))
end


function int_nlsolve!(dx,x,zp)
    z, p = zp[1:M], zp[M+1:6N+1+2N*M+M^2]
    c, α, β = selc*z, selα*z, selβ*z
    obj_x(x) = obj_nlsolve(z,p,x)
    dx .= (obj_x(x).*extpow_prod_pdf(x,α,β))[1,:]
end


function int_cont!(dx,x,zp)
  z, p = zp[1:M], zp[M+1:6N+1+2N*M+M^2]
  c, α, β = selc*z, selα*z, selβ*z
  obj_x(x) = obj_contraction(z,p,x)
  dx .= (obj_x(x).* extpow_prod_pdf(x,α,β))[1,:]
end

function int_nlsolve(x,zp)
    z, p = zp[1:M], zp[M+1:6N+1+2N*M+M^2]
    c, α, β = selc*z, selα*z, selβ*z
    obj_x(x) = obj_nlsolve(z,p,x)
    (obj_x(x).* extpow_prod_pdf(x,α,β))[1,:]
end

function int_cont(x,zp)
  z, p = zp[1:M], zp[M+1:6N+1+2N*M+M^2]
  c, α, β = selc*z, selα*z, selβ*z
  obj_x(x) = obj_contraction(z,p,x)
  (obj_x(x).*extpow_prod_pdf(x,α,β))[1,:]
end

# function obj_nlsolve_int(z,p)  
#     mb = size(selp,1)÷N
#     prob = QuadratureProblem(int_nlsolve!,zeros(N),ones(N),vcat(z,p),batch=mb)
#     sol = Quadrature.solve(prob,CubaCuhre(),reltol=1e-2,abstol=1e-2)[1]
# end

# function obj_cont_int(z,p)
#     mb = size(selp,1)÷N
#     prob = QuadratureProblem(int_cont!,zeros(N),ones(N),vcat(z,p),batch=mb)
#     sol = Quadrature.solve(prob,CubaCuhre(),reltol=1e-1,abstol=1e-1)[1]
# end

# function obj_nlsolve_intH(z,p)
#     mb = size(selp,1)÷N
#     prob = QuadratureProblem(int_nlsolve,zeros(N),ones(N),vcat(z,p))
#     sol = Quadrature.solve(prob,HCubatureJL(),reltol=1e-1,abstol=1e-1)[1] # CubatureJLp, CubatureJLh, HCubatureJL
# end

# function obj_cont_intH(z,p)
#     mb = size(selp,1)÷N
#     prob = QuadratureProblem(int_cont,zeros(N),ones(N),vcat(z,p))
#     sol = Quadrature.solve(prob,HCubatureJL(),reltol=1e-1,abstol=1e-1)[1] # CubatureJLp, CubatureJLh, HCubatureJL
# end

# upper and lower bounds
# low<=z<=upp
lowα = 0.01 ;uppα = 0.5
lowβ = 0.01   ;uppβ = 0.5
low_zero = eps(T) # eps(T) T(0.005)
low = vcat(spzeros(T,N^2),fill(low_zero,N),fill(low_zero,N),fill(lowα,N),fill(lowβ,N))
lowA = Array(low)
upp = vcat(ones(T,N^2),p_bar,a,Array{T}(fill(uppα,N)),Array{T}(fill(uppβ,N)))

const p_sparse = vcat(sparse(vcat(w, p_bar, delta, γ, a ,q)),reshape(hᵀ,2N*M),reshape(H,M^2))
#const p_dense =  vcat(w, p_bar, delta, γ, a ,q,Array(hᵀ)...,Array(H)...)
P = 2N*M+M^2+6N+1 #P==length(p_dense)==number of parameters

z0 = sparse(lowA + (upp-lowA)/2);
const A0 = (LowerTriangular(ones(N,N))-Matrix(I,N,N))/2
z0[1:N^2] = [A0...]
const p0 = p_sparse #p_dense

#x0 = ones(N,10)/10 #ones(N,mini_batch)/10
# rescale_x(x) = (x .+ 1)./2
# quadInt(f) = dot(map(x -> f(rescale_x(x)), nodes0), weights0)/2^N
# nodes1 = map(x->(x .+ 1) ./ 2,nodes0)
# quadInt2(f) = dot(map(x -> f(x), nodes1), weights0)/2^N
#@test quadInt(x->(sum(x)-2.5)^2) ≈ 5/12

end # module